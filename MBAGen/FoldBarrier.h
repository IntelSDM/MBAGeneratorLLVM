#pragma once

#include "llvm/IR/BasicBlock.h"
#include "llvm/IR/Constants.h"
#include "llvm/IR/DerivedTypes.h"
#include "llvm/IR/Function.h"
#include "llvm/IR/GlobalVariable.h"
#include "llvm/IR/IRBuilder.h"
#include "llvm/IR/InlineAsm.h"
#include "llvm/IR/Instructions.h"
#include "llvm/IR/Metadata.h"
#include "llvm/IR/Module.h"
#include "llvm/IR/Type.h"
#include "llvm/IR/Value.h"
#include "llvm/Transforms/Utils/ModuleUtils.h"

#include <cassert>
#include <cstdint>
#include <random>
#include <string>
#include <vector>

namespace mbagen
{

enum class BarrierShape : uint8_t
{
    Legacy = 0,
    XorHidden = 1,
    AddHidden = 2,
    SubHidden = 3,
};

struct BarrierSlot
{
    llvm::GlobalVariable* GV = nullptr;
    llvm::GlobalVariable* NonceGV = nullptr;
    llvm::GlobalVariable* CipherGV = nullptr;
    uint64_t K = 0;
    uint64_t Nonce = 0;
    uint64_t Cipher = 0;
};

struct ModuleBarrierState
{
    llvm::Module* M = nullptr;
    unsigned PoolSize = 16;
    std::vector<BarrierSlot> SlotsI32;
    std::vector<BarrierSlot> SlotsI64;
    std::mt19937_64 Rng;
    unsigned WLegacy = 0;
    unsigned WXor = 50;
    unsigned WAdd = 25;
    unsigned WSub = 25;
    std::string SlotNamePrefix = ".mba.barrier";

    llvm::GlobalVariable* SeedGV = nullptr;
    llvm::GlobalVariable* SeedBootNonceGV = nullptr;
    llvm::GlobalVariable* SeedBootSaltGV = nullptr;
    uint64_t SeedK = 0;
    uint64_t SeedBootNonce = 0;
    uint64_t SeedBootSalt = 0;
    bool InitEmitted = false;

    void Init(llvm::Module* module, unsigned poolsz, uint64_t seed, unsigned wleg, unsigned wx, unsigned wa, unsigned ws, const std::string& nameprefix)
    {
        M = module;
        PoolSize = (poolsz == 0) ? 1u : poolsz;
        SlotsI32.assign(PoolSize, BarrierSlot{});
        SlotsI64.assign(PoolSize, BarrierSlot{});
        // seed=0 is reserved for "no seed given, derive one"; using the module pointer XOR'd with the golden ratio constant keeps builds deterministic within a translation unit without colliding across modules whose pointers happen to agree in some bits
        Rng.seed(seed ? seed : (uint64_t)(uintptr_t)module ^ 0x9E3779B97F4A7C15ull);
        WLegacy = wleg;
        WXor = wx;
        WAdd = wa;
        WSub = ws;
        SlotNamePrefix = nameprefix;

        SeedGV = nullptr;
        SeedBootNonceGV = nullptr;
        SeedBootSaltGV = nullptr;
        // both bootstrap halves forced odd so the first ARX round's shift by 13 feeds a non degenerate low bit; a zero bootstrap would collapse the schedule to (0, 0) across every slot and leak K == Cipher as a plain XOR with a fixed constant
        SeedBootNonce = Rng() | 1ull;
        SeedBootSalt = Rng() | 1ull;
        SeedK = 0;
        InitEmitted = false;
    }
};

// ARX schedule must round trip bit-identical between pass time Cipher precompute and the runtime decrypt emitted by BarrierEmitArxBody; any divergence between this C++ and that IR silently poisons every K slot and the whole barrier pool decodes to garbage on first launch
static inline uint64_t BarrierArxMix(uint64_t nonce, uint64_t seed)
{
    uint64_t x = nonce ^ seed;
    x = ((x << 13) | (x >> 51)) + 0x9E3779B97F4A7C15ULL;
    x = x ^ ((x << 7) | (x >> 57));
    x = ((x << 17) | (x >> 47)) + seed;
    x = x ^ nonce;
    x = ((x << 23) | (x >> 41)) * 0xBF58476D1CE4E5B9ULL;
    return x;
}

static BarrierShape BarrierPickShape(ModuleBarrierState& st)
{
    const unsigned total = st.WLegacy + st.WXor + st.WAdd + st.WSub;
    if (total == 0) return BarrierShape::Legacy;
    unsigned r = (unsigned)(st.Rng() % total);
    if (r < st.WLegacy) return BarrierShape::Legacy;
    r -= st.WLegacy;
    if (r < st.WXor) return BarrierShape::XorHidden;
    r -= st.WXor;
    if (r < st.WAdd) return BarrierShape::AddHidden;
    return BarrierShape::SubHidden;
}

static uint64_t BarrierDrawK(ModuleBarrierState& st, unsigned w)
{
    const uint64_t mask = (w >= 64) ? ~uint64_t(0) : ((uint64_t(1) << w) - 1);
    uint64_t k = 0;
    // avoid k in {0, 1, 2} and their all ones duals; those shapes decompile as near identity (xor with zero, add with one) and give an attacker eyeballing the disassembly a free hint about which operations are barriers vs real work
    do
    {
        k = st.Rng() & mask;
    } while (k == 0 || k == 1 || k == 2 || k == mask || k == (mask - 1) || k == (mask - 2));
    return k;
}

static llvm::GlobalVariable* BarrierCreateConstDataGV(llvm::Module* m, llvm::IntegerType* ity, uint64_t value, const std::string& name)
{
    using namespace llvm;
    // const data operand GVs are NOT constant at the LLVM level; marking them constant lets GVN propagate the stored initializer into every volatile load inside the init function and the ARX computation folds at compile time, collapsing the whole scheme back to a literal K
    GlobalVariable* gv = new GlobalVariable(*m, ity, /*isConstant=*/false, GlobalValue::InternalLinkage, ConstantInt::get(ity, value), name);
    gv->setUnnamedAddr(GlobalValue::UnnamedAddr::None);
    gv->setAlignment(Align(ity->getBitWidth() / 8));
    return gv;
}

static void BarrierEnsureSeedGVs(ModuleBarrierState& st)
{
    using namespace llvm;
    if (st.SeedGV) return;
    LLVMContext& ctx = st.M->getContext();
    IntegerType* i64 = IntegerType::get(ctx, 64);

    st.SeedK = BarrierArxMix(st.SeedBootNonce, st.SeedBootSalt);

    // SeedGV is zero initialised so the on disk .bss (not .data) image reveals nothing; the runtime initialiser recomputes SeedK via BarrierArxMix(SeedBootNonce, SeedBootSalt) and stores it here before any protected function runs
    st.SeedGV = new GlobalVariable(*st.M, i64, /*isConstant=*/false, GlobalValue::InternalLinkage, ConstantInt::get(i64, 0), st.SlotNamePrefix + ".seed");
    st.SeedGV->setUnnamedAddr(GlobalValue::UnnamedAddr::None);
    st.SeedGV->setAlignment(Align(8));

    st.SeedBootNonceGV = BarrierCreateConstDataGV(st.M, i64, st.SeedBootNonce, st.SlotNamePrefix + ".seed.nonce");
    st.SeedBootSaltGV = BarrierCreateConstDataGV(st.M, i64, st.SeedBootSalt, st.SlotNamePrefix + ".seed.salt");
}

static BarrierSlot& BarrierGetOrCreateSlot(ModuleBarrierState& st, llvm::IntegerType* ity, unsigned idx)
{
    using namespace llvm;
    const unsigned w = ity->getBitWidth();
    auto& pool = (w == 32) ? st.SlotsI32 : st.SlotsI64;
    unsigned i = idx % (unsigned)pool.size();
    BarrierSlot& slot = pool[i];
    if (slot.GV) return slot;

    BarrierEnsureSeedGVs(st);

    slot.K = BarrierDrawK(st, w);
    // Nonce forced odd so ARX(Nonce, Seed) across slots cannot collide via a zero nonce shared with another path; the entropy loss is one bit and the uniqueness guarantee is per pool
    slot.Nonce = st.Rng() | 1ull;
    const uint64_t derived = BarrierArxMix(slot.Nonce, st.SeedK);
    const uint64_t wmask = (w >= 64) ? ~uint64_t(0) : ((uint64_t(1) << w) - 1);
    slot.Cipher = (slot.K ^ derived) & wmask;

    const std::string base = st.SlotNamePrefix + ".g" + std::to_string(w) + "." + std::to_string(i);

    slot.GV = new GlobalVariable(*st.M, ity, /*isConstant=*/false, GlobalValue::InternalLinkage, ConstantInt::get(ity, 0), base);
    // setUnnamedAddr(None) is required; UnnamedAddr::Global would merge two independent slots whose K draws happened to collide, collapsing the pool and reusing the same GV for sites that should see independent hidden constants
    slot.GV->setUnnamedAddr(GlobalValue::UnnamedAddr::None);
    slot.GV->setAlignment(Align(w / 8));

    IntegerType* i64 = IntegerType::get(st.M->getContext(), 64);
    slot.NonceGV = BarrierCreateConstDataGV(st.M, i64, slot.Nonce, base + ".nonce");
    slot.CipherGV = BarrierCreateConstDataGV(st.M, ity, slot.Cipher, base + ".cipher");

    return slot;
}

static unsigned BarrierPickSecondSlotIndex(ModuleBarrierState& st, unsigned firstidx)
{
    const unsigned pool = st.PoolSize;
    if (pool < 2) return firstidx;
    // offset is in [1, pool-1] so (firstidx + off) % pool never lands back on firstidx; sampling (rng() % pool) and retrying on equality would skew the distribution toward the hole in pathological pool sizes
    const unsigned off = 1u + (unsigned)(st.Rng() % (pool - 1u));
    return (firstidx + off) % pool;
}

static void BarrierTagConstBarrierEscalate(llvm::IRBuilder<>& builder, llvm::Value* v)
{
    using namespace llvm;
    auto* inst = dyn_cast<Instruction>(v);
    if (!inst) return;
    LLVMContext& ctx = builder.getContext();
    MDNode* md = MDNode::get(ctx, { MDString::get(ctx, "const_barrier") });
    // the "mbagen.escalate" metadata key paired with "const_barrier" is the downstream contract ConstMBAPass reads to route k at the barrier site into its escalated tier MBA emission; renaming either side silently breaks the compounding defense
    inst->setMetadata("mbagen.escalate", md);
}

static llvm::LoadInst* BarrierEmitVolatileLoad(llvm::IRBuilder<>& builder, llvm::GlobalVariable* gv, llvm::IntegerType* ity, const char* name)
{
    using namespace llvm;
    // load MUST be volatile; a non volatile load lets LLVM fold through the stored initial value of the GV (zero for K slots, the plaintext constant for Nonce/Cipher GVs), collapsing the barrier at compile time and defeating the whole hidden constant construction
    LoadInst* ld = builder.CreateLoad(ity, gv, /*isVolatile=*/true, name);
    ld->setAlignment(Align(ity->getBitWidth() / 8));
    return ld;
}

inline llvm::Value* InsertFoldBarrierLegacy(llvm::IRBuilder<>& builder, llvm::Value* v)
{
    using namespace llvm;
    Type* ty = v->getType();
    assert(ty->isIntegerTy() && "FoldBarrier only supports integer types");
    FunctionType* ft = FunctionType::get(ty, { ty }, /*isVarArg=*/false);
    InlineAsm* ia = InlineAsm::get(ft, /*AsmString=*/"", /*Constraints=*/"=r,0", /*hasSideEffects=*/true, /*isAlignStack=*/false, InlineAsm::AD_ATT);
    return builder.CreateCall(ia, { v }, "MBABarrierLeg");
}

inline llvm::Value* InsertFoldBarrier(llvm::IRBuilder<>& builder, llvm::Value* v, ModuleBarrierState& st)
{
    using namespace llvm;
    Type* ty = v->getType();
    assert(ty->isIntegerTy() && "FoldBarrier only supports integer types");
    auto* ity = cast<IntegerType>(ty);
    const unsigned w = ity->getBitWidth();

    // slot pools are maintained only for i32 and i64 because ConstMBA/LinearMBA/PolyMBA all emit at these widths; any other width falls cleanly through to the legacy inline asm shape rather than materialising a half configured pool at that width
    if (w != 32 && w != 64) return InsertFoldBarrierLegacy(builder, v);

    BarrierShape shape = BarrierPickShape(st);
    if (shape == BarrierShape::Legacy) return InsertFoldBarrierLegacy(builder, v);

    const unsigned slotidxa = (unsigned)(st.Rng() % st.PoolSize);
    const unsigned slotidxb = BarrierPickSecondSlotIndex(st, slotidxa);
    BarrierSlot& slota = BarrierGetOrCreateSlot(st, ity, slotidxa);
    BarrierSlot& slotb = (st.PoolSize >= 2 && slotidxb != slotidxa)
        ? BarrierGetOrCreateSlot(st, ity, slotidxb)
        : slota;

    const bool dualslot = (&slota != &slotb);
    const uint64_t wmask = (w >= 64) ? ~uint64_t(0) : ((uint64_t(1) << w) - 1);

    LoadInst* lda = BarrierEmitVolatileLoad(builder, slota.GV, ity, "MBABarrierLdA");
    LoadInst* ldb = dualslot ? BarrierEmitVolatileLoad(builder, slotb.GV, ity, "MBABarrierLdB") : nullptr;

    Value* result = nullptr;
    switch (shape)
    {
        case BarrierShape::XorHidden:
        {
            Value* t1 = builder.CreateXor(v, lda, "MBABarrierX1");
            if (dualslot) t1 = builder.CreateXor(t1, ldb, "MBABarrierX2");
            
            const uint64_t kadj = dualslot ? ((slota.K ^ slotb.K) & wmask) : (slota.K & wmask);
            Value* kc = ConstantInt::get(ity, kadj);
            Value* t2 = builder.CreateXor(t1, kc, "MBABarrierXor");
            BarrierTagConstBarrierEscalate(builder, t2);
            result = t2;
            break;
        }
        case BarrierShape::AddHidden:
        {
            Value* t1 = builder.CreateAdd(v, lda, "MBABarrierA1");
            if (dualslot) t1 = builder.CreateAdd(t1, ldb, "MBABarrierA2");
            const uint64_t kadj = dualslot ? ((slota.K + slotb.K) & wmask) : (slota.K & wmask);
            Value* kc = ConstantInt::get(ity, kadj);
            Value* t2 = builder.CreateSub(t1, kc, "MBABarrierAdd");
            BarrierTagConstBarrierEscalate(builder, t2);
            result = t2;
            break;
        }
        case BarrierShape::SubHidden:
        {
            Value* t1 = builder.CreateSub(v, lda, "MBABarrierS1");
            if (dualslot) t1 = builder.CreateSub(t1, ldb, "MBABarrierS2");
            const uint64_t kadj = dualslot ? ((slota.K + slotb.K) & wmask) : (slota.K & wmask);
            Value* kc = ConstantInt::get(ity, kadj);
            Value* t2 = builder.CreateAdd(t1, kc, "MBABarrierSub");
            BarrierTagConstBarrierEscalate(builder, t2);
            result = t2;
            break;
        }
        case BarrierShape::Legacy:
        default:
            return InsertFoldBarrierLegacy(builder, v);
    }

    return result;
}

inline llvm::Value* InsertFoldBarrier(llvm::IRBuilder<>& builder, llvm::Value* v)
{
    return InsertFoldBarrierLegacy(builder, v);
}

// step sequence here must stay 1:1 with BarrierArxMix; a swap or an omitted step produces a different cipher to K relation and the whole module decrypts to zero on first launch
static llvm::Value* BarrierEmitArxBody(llvm::IRBuilder<>& builder, llvm::Value* nonce, llvm::Value* seed)
{
    using namespace llvm;
    LLVMContext& ctx = builder.getContext();
    IntegerType* i64 = IntegerType::get(ctx, 64);

    auto rotl = [&](Value* v, unsigned n) -> Value*
    {
        Value* lo = builder.CreateShl(v, ConstantInt::get(i64, n));
        Value* hi = builder.CreateLShr(v, ConstantInt::get(i64, 64 - n));
        return builder.CreateOr(lo, hi);
    };

    Value* x = builder.CreateXor(nonce, seed);
    x = builder.CreateAdd(rotl(x, 13), ConstantInt::get(i64, 0x9E3779B97F4A7C15ULL));
    x = builder.CreateXor(x, rotl(x, 7));
    x = builder.CreateAdd(rotl(x, 17), seed);
    x = builder.CreateXor(x, nonce);
    x = builder.CreateMul(rotl(x, 23), ConstantInt::get(i64, 0xBF58476D1CE4E5B9ULL));
    return x;
}

static void BarrierEmitInitBody(llvm::IRBuilder<>& builder, ModuleBarrierState& st)
{
    using namespace llvm;
    LLVMContext& ctx = builder.getContext();
    IntegerType* i64 = IntegerType::get(ctx, 64);

    // seed derivation runs first and the seed store is reloaded via a fresh volatile load before any slot store reads it; without that reload the IRBuilder could reuse the pre store SSA value and skip the memory round trip that forces LLVM to treat Seed as runtime defined
    LoadInst* bootnoncev = BarrierEmitVolatileLoad(builder, st.SeedBootNonceGV, i64, "wp.bar.sbn");
    LoadInst* bootsaltv = BarrierEmitVolatileLoad(builder, st.SeedBootSaltGV, i64, "wp.bar.sbs");
    Value* seedval = BarrierEmitArxBody(builder, bootnoncev, bootsaltv);
    StoreInst* seedstore = builder.CreateStore(seedval, st.SeedGV);
    seedstore->setVolatile(true);
    seedstore->setAlignment(Align(8));

    LoadInst* seedld = BarrierEmitVolatileLoad(builder, st.SeedGV, i64, "wp.bar.seed");

    IntegerType* i32 = IntegerType::get(ctx, 32);

    auto emitslotstore = [&](const BarrierSlot& slot, IntegerType* ity)
    {
        if (!slot.GV || !slot.NonceGV || !slot.CipherGV) return;
        const unsigned w = ity->getBitWidth();
        LoadInst* nonceld = BarrierEmitVolatileLoad(builder, slot.NonceGV, i64, "wp.bar.n");
        Value* derived = BarrierEmitArxBody(builder, nonceld, seedld);
        if (w != 64) derived = builder.CreateTrunc(derived, ity, "wp.bar.dtr");
        LoadInst* cipherld = BarrierEmitVolatileLoad(builder, slot.CipherGV, ity, "wp.bar.c");
        Value* kval = builder.CreateXor(cipherld, derived, "wp.bar.k");
        StoreInst* kst = builder.CreateStore(kval, slot.GV);
        kst->setVolatile(true);
        kst->setAlignment(Align(w / 8));
    };

    for (const BarrierSlot& s : st.SlotsI64) emitslotstore(s, i64);
    for (const BarrierSlot& s : st.SlotsI32) emitslotstore(s, i32);
}

inline void BarrierFinalizeModule(ModuleBarrierState& st)
{
    using namespace llvm;
    if (!st.M || st.InitEmitted) return;
    bool anyslot = false;
    for (const auto& s : st.SlotsI32) if (s.GV) { anyslot = true; break; }
    if (!anyslot) for (const auto& s : st.SlotsI64) if (s.GV) { anyslot = true; break; }
    // no slots touched => no seed GV was created either, so emitting an init that stores to a null SeedGV would crash; the no op early out keeps RunConstMBA/RunLinearMBA/RunPolyMBA callable on modules that never materialised a hidden shape barrier
    if (!anyslot) return;

    LLVMContext& ctx = st.M->getContext();
    FunctionType* ft = FunctionType::get(Type::getVoidTy(ctx), /*isVarArg=*/false);
    const std::string fname = "__wp_barrier_init" + st.SlotNamePrefix;
    Function* initfn = Function::Create(ft, GlobalValue::InternalLinkage, fname, st.M);
    // bug: NoInline keeps the init body out of llvm.global_ctors' thunk and out of any caller that happens to pick it up through LTO; OptimizeNone was previously set as belt-and-braces against GVN propagating stored initializers through the volatile loads but caused downstream pipeline crashes on PolyMBA-heavy modules, so it's dropped — the volatile load/store semantics are already sufficient to block the fold
    initfn->addFnAttr(Attribute::NoInline);
    BasicBlock* entry = BasicBlock::Create(ctx, "Entry", initfn);
    IRBuilder<> ib(entry);
    BarrierEmitInitBody(ib, st);
    ib.CreateRetVoid();

    appendToGlobalCtors(*st.M, initfn, 101);

    st.InitEmitted = true;
}

}