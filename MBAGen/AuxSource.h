#pragma once

#include "llvm/ADT/DenseMap.h"
#include "llvm/ADT/SmallVector.h"
#include "llvm/IR/Argument.h"
#include "llvm/IR/BasicBlock.h"
#include "llvm/IR/Constants.h"
#include "llvm/IR/Dominators.h"
#include "llvm/IR/Function.h"
#include "llvm/IR/IRBuilder.h"
#include "llvm/IR/Instructions.h"
#include "llvm/IR/Type.h"
#include "llvm/IR/Value.h"

#include <algorithm>
#include <random>
#include <vector>

namespace mbagen
{

inline constexpr unsigned AuxPoolSlots = 8;

inline llvm::DenseMap<llvm::Function*, llvm::AllocaInst*> gAuxPoolCache;

inline void ResetAuxPoolCache()
{
    gAuxPoolCache.clear();
}

inline llvm::AllocaInst* GetOrCreateAuxPool(llvm::Function& f, std::mt19937_64& rng)
{
    using namespace llvm;

    auto it = gAuxPoolCache.find(&f);
    if (it != gAuxPoolCache.end()) return it->second;

    BasicBlock& entry = f.getEntryBlock();
    BasicBlock::iterator insertpt = entry.getFirstInsertionPt();

    LLVMContext& ctx = f.getContext();
    IRBuilder<> builder(&entry, insertpt);

    Type* i64 = Type::getInt64Ty(ctx);
    ArrayType* poolty = ArrayType::get(i64, AuxPoolSlots);

    AllocaInst* pool = builder.CreateAlloca(poolty, nullptr, ".mbagen.aux_pool");
    pool->setAlignment(Align(8));

    Type* i32 = Type::getInt32Ty(ctx);
    for (unsigned i = 0; i < AuxPoolSlots; ++i)
    {
        Value* idx[2] = { ConstantInt::get(i32, 0), ConstantInt::get(i32, i) };
        Value* slotptr = builder.CreateInBoundsGEP(poolty, pool, idx, ".mbagen.aux_slot");
        uint64_t init = rng();
        // stores MUST be volatile; without the volatile flag, InstCombine propagates the stored initial constants through to every volatile load at use sites, collapsing the aux pool at compile time and leaving the MBA expression as a visible constant
        StoreInst* st = builder.CreateStore(ConstantInt::get(i64, init), slotptr);
        st->setVolatile(true);
        st->setAlignment(Align(8));
    }

    gAuxPoolCache[&f] = pool;
    return pool;
}

inline std::vector<llvm::Value*> HarvestLiveIntegerValues(llvm::Instruction* insertionpoint, const llvm::DominatorTree& dt, llvm::IntegerType* requiredtype, unsigned maxcount)
{
    using namespace llvm;

    std::vector<Value*> candidates;
    candidates.reserve(maxcount);

    Function* f = insertionpoint->getFunction();
    if (!f) return candidates;

    for (Argument& a : f->args())
    {
        if (a.getType() != requiredtype) continue;
        candidates.push_back(&a);
        if (candidates.size() >= maxcount) return candidates;
    }

    for (BasicBlock& bb : *f)
    {
        for (Instruction& i : bb)
        {
            if (&i == insertionpoint) continue;
            if (i.getType() != requiredtype) continue;
            // domination check rejects not just self reference but also PHIs with incomplete incoming and the insertion point itself; without it a replace operand use would self reference and the IR would silently build a cycle
            if (!dt.dominates(&i, insertionpoint)) continue;
            candidates.push_back(&i);
            if (candidates.size() >= maxcount) return candidates;
        }
    }

    return candidates;
}

inline llvm::Value* EmitAuxPoolLoad(llvm::IRBuilder<>& builder, llvm::AllocaInst* pool, llvm::IntegerType* resulttype, unsigned slotidx)
{
    using namespace llvm;

    LLVMContext& ctx = builder.getContext();
    Type* i32 = Type::getInt32Ty(ctx);
    Type* i64 = Type::getInt64Ty(ctx);

    ArrayType* poolty = ArrayType::get(i64, AuxPoolSlots);
    Value* idx[2] = { ConstantInt::get(i32, 0), ConstantInt::get(i32, slotidx % AuxPoolSlots) };
    Value* slotptr = builder.CreateInBoundsGEP(poolty, pool, idx, ".mbagen.aux_slot");

    LoadInst* ld = builder.CreateLoad(i64, slotptr, ".mbagen.aux");
    ld->setVolatile(true);
    ld->setAlignment(Align(8));

    if (resulttype == i64) return ld;
    // the trunc here sits in the backfill path only, it does NOT contaminate the MBA shape proper (which sees a plain integer-typed Value*); auto zext/sext/trunc on the Tier A harvest path was deliberately excluded for exactly this reason
    return builder.CreateTrunc(ld, resulttype, ".mbagen.aux_tr");
}

inline std::vector<llvm::Value*> GatherAuxValues(llvm::IRBuilder<>& builder, const llvm::DominatorTree& dt, llvm::IntegerType* resulttype, unsigned V, std::mt19937_64& rng, unsigned crosswidthpct = 0)
{
    using namespace llvm;

    if (crosswidthpct > 100) crosswidthpct = 100;

    Instruction* insertpt = &*builder.GetInsertPoint();
    std::vector<Value*> harvested = HarvestLiveIntegerValues(insertpt, dt, resulttype, V * 3);

    std::vector<Value*> crossharvested;
    IntegerType* crosstype = nullptr;
    if (crosswidthpct > 0)
    {
        LLVMContext& ctx = builder.getContext();
        if (resulttype->getBitWidth() == 32) crosstype = Type::getInt64Ty(ctx);
        else if (resulttype->getBitWidth() == 64) crosstype = Type::getInt32Ty(ctx);
        if (crosstype)
        {
            crossharvested = HarvestLiveIntegerValues(insertpt, dt, crosstype, V * 3);
            if (crossharvested.size() > 1) std::shuffle(crossharvested.begin(), crossharvested.end(), rng);
        }
    }

    // function arguments always rank first in the harvest (dominance by construction), and without the shuffle that ordering becomes a static tell "first V aux slots are always function args" so an attacker can scan for exactly that header in each MBA site
    if (harvested.size() > 1) std::shuffle(harvested.begin(), harvested.end(), rng);

    std::vector<Value*> result;
    result.reserve(V);

    const unsigned fromharvest = std::min(V, (unsigned)harvested.size());
    size_t crosscursor = 0;
    unsigned samecursor = 0;
    while (result.size() < fromharvest || (crosswidthpct > 0 && crosscursor < crossharvested.size() && result.size() < V))
    {
        if (result.size() >= V) break;
        const bool havesame = (samecursor < harvested.size());
        const bool havecross = (crosswidthpct > 0 && crosscursor < crossharvested.size());
        const bool rollcross = havecross && ((unsigned)(rng() % 100u) < crosswidthpct);
        if (rollcross)
        {
            Value* cv = crossharvested[crosscursor++];
            Value* bridged;
            if (crosstype->getBitWidth() > resulttype->getBitWidth())
                bridged = builder.CreateTrunc(cv, resulttype, ".mbagen.aux_xw_tr");
            else
                bridged = builder.CreateZExt(cv, resulttype, ".mbagen.aux_xw_zx");
            result.push_back(bridged);
        }
        else if (havesame)
        {
            result.push_back(harvested[samecursor++]);
        }
        else
        {
            break;
        }
    }

    if (result.size() < V)
    {
        Function* f = insertpt->getFunction();
        AllocaInst* pool = GetOrCreateAuxPool(*f, rng);

        unsigned slot = (unsigned)(rng() % AuxPoolSlots);
        while (result.size() < V)
        {
            result.push_back(EmitAuxPoolLoad(builder, pool, resulttype, slot));
            slot = (slot + 1) % AuxPoolSlots;
        }
    }

    return result;
}

}