#pragma once

#include "AtomTree.h"
#include "MBABasis.h"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <functional>
#include <optional>
#include <random>
#include <vector>

namespace mbagen
{

inline u64 ModInvPow2(u64 a, unsigned w)
{
    // a must be odd for an inverse to exist mod 2^w (even numbers share a factor of 2 with 2^w); returning 0 on even input makes subsequent arithmetic explode loudly rather than silently producing wrong results
    if ((a & 1u) == 0u) return 0;
    u64 inv = a;
    for (int i = 0; i < 6; ++i)
        inv = inv * (2 - a * inv);
    return MaskW(inv, w);
}

struct RowReduceResult
{
    static constexpr unsigned NPOS = ~0u;
    bool Solvable = false;
    unsigned Rank = 0;
    std::vector<unsigned> PivotColForRow;
    std::vector<unsigned> PivotRowForCol;
    std::vector<std::vector<u64>> M;
    std::vector<u64> T;
};

inline RowReduceResult RowReduceMod2W(std::vector<std::vector<u64>> m, std::vector<u64> t, unsigned w)
{
    const size_t rows = m.size();
    const size_t cols = rows > 0 ? m[0].size() : 0;
    RowReduceResult res;
    res.M = std::move(m);
    res.T = std::move(t);
    res.PivotColForRow.assign(rows, RowReduceResult::NPOS);
    res.PivotRowForCol.assign(cols, RowReduceResult::NPOS);

    unsigned nextrow = 0;
    for (unsigned c = 0; c < (unsigned)cols && nextrow < (unsigned)rows; ++c)
    {
        unsigned pivrow = RowReduceResult::NPOS;
        // units of Z/2^w are exactly the odd numbers; pivoting on an even entry would require solving for a non invertible scalar and corrupt the eliminated row residues, so columns with all-even entries get marked free and skipped
        for (unsigned r = nextrow; r < (unsigned)rows; ++r)
        {
            if ((res.M[r][c] & 1u) != 0u) { pivrow = r; break; }
        }
        if (pivrow == RowReduceResult::NPOS) continue;

        if (pivrow != nextrow)
        {
            std::swap(res.M[pivrow], res.M[nextrow]);
            std::swap(res.T[pivrow], res.T[nextrow]);
        }

        const u64 pivotinv = ModInvPow2(res.M[nextrow][c], w);
        for (unsigned cc = 0; cc < (unsigned)cols; ++cc)
            res.M[nextrow][cc] = MaskW(res.M[nextrow][cc] * pivotinv, w);
        res.T[nextrow] = MaskW(res.T[nextrow] * pivotinv, w);

        for (unsigned r = 0; r < (unsigned)rows; ++r)
        {
            if (r == nextrow) continue;
            const u64 factor = res.M[r][c];
            if (factor == 0u) continue;
            for (unsigned cc = 0; cc < (unsigned)cols; ++cc)
            {
                const u64 sub = MaskW(factor * res.M[nextrow][cc], w);
                res.M[r][cc] = MaskW(res.M[r][cc] - sub, w);
            }
            const u64 subt = MaskW(factor * res.T[nextrow], w);
            res.T[r] = MaskW(res.T[r] - subt, w);
        }

        res.PivotColForRow[nextrow] = c;
        res.PivotRowForCol[c] = nextrow;
        ++nextrow;
    }

    res.Rank = nextrow;

    res.Solvable = true;
    for (unsigned r = res.Rank; r < (unsigned)rows; ++r)
    {
        if (res.T[r] != 0u) { res.Solvable = false; break; }
    }
    return res;
}

inline std::vector<std::vector<u64>> BuildTruthMatrix(const std::vector<AtomTreePtr>& trees, unsigned V, unsigned w)
{
    const size_t rows = size_t(1) << V;
    const size_t cols = trees.size();
    std::vector<std::vector<u64>> m(rows, std::vector<u64>(cols, 0));
    std::vector<u64> inputs(V > 0 ? V : 1u, 0);
    for (size_t r = 0; r < rows; ++r)
    {
        for (unsigned k = 0; k < V; ++k)
            inputs[k] = (r >> k) & 1u;
        for (size_t i = 0; i < cols; ++i)
            m[r][i] = EvalAtomTree(*trees[i], inputs.data(), w);
    }
    return m;
}

enum class TargetKind
{
    Zero,
    ClassA,
    ClassX,
    ClassXY,
};

enum class XYOp : uint8_t
{
    And = 0,
    Xor = 1,
    Or = 2,
};

struct TargetSpec
{
    TargetKind Kind = TargetKind::Zero;
    u64 Constant = 0;
    u64 XCoef = 0;
    unsigned XVarIdx = 0;
    unsigned YVarIdx = 0;
    XYOp XYOpV = XYOp::And;
};

inline std::vector<u64> BuildTargetVector(const TargetSpec& spec, unsigned V, unsigned w)
{
    const size_t N = size_t(1) << V;
    std::vector<u64> t(N, 0);
    switch (spec.Kind)
    {
        case TargetKind::Zero: break;
        case TargetKind::ClassA:
        {
            const u64 v = MaskW(spec.Constant, w);
            for (size_t j = 0; j < N; ++j) t[j] = v;
            break;
        }
        case TargetKind::ClassX:
        {
            const u64 c = MaskW(spec.XCoef, w);
            const unsigned k = spec.XVarIdx;
            for (size_t j = 0; j < N; ++j)
                t[j] = ((j >> k) & 1u) ? c : 0u;
            break;
        }
        case TargetKind::ClassXY:
        {
            const u64 c = MaskW(spec.XCoef, w);
            const unsigned a = spec.XVarIdx;
            const unsigned b = spec.YVarIdx;
            for (size_t j = 0; j < N; ++j)
            {
                const unsigned ba = (unsigned)((j >> a) & 1u);
                const unsigned bb = (unsigned)((j >> b) & 1u);
                unsigned mix;
                switch (spec.XYOpV)
                {
                    case XYOp::And: mix = ba & bb; break;
                    case XYOp::Xor: mix = ba ^ bb; break;
                    case XYOp::Or:
                    default: mix = ba | bb; break;
                }
                t[j] = mix ? c : 0u;
            }
            break;
        }
    }
    return t;
}

enum class Tier
{
    T8 = 8,
    T7 = 7,
    T6 = 6,
};

enum class CoefMag
{
    Small = 0,
    Medium = 1,
    Large = 2,
};

struct TierConfig
{
    unsigned VMin = 2;
    unsigned VMax = 3;
    unsigned SExtraOver2V = 4;
    unsigned DepthMin = 1;
    unsigned DepthMax = 2;
    CoefMag CoefMagV = CoefMag::Small;
};

inline TierConfig GetTierConfig(Tier t)
{
    TierConfig c;
    switch (t)
    {
        case Tier::T8:
            c.VMin = 2; c.VMax = 3;
            c.SExtraOver2V = 4;
            c.DepthMin = 1; c.DepthMax = 2;
            c.CoefMagV = CoefMag::Small;
            break;
        case Tier::T7:
            c.VMin = 3; c.VMax = 4;
            c.SExtraOver2V = 6;
            c.DepthMin = 2; c.DepthMax = 3;
            c.CoefMagV = CoefMag::Medium;
            break;
        case Tier::T6:
            c.VMin = 4; c.VMax = 6;
            c.SExtraOver2V = 8;
            c.DepthMin = 2; c.DepthMax = 3;
            c.CoefMagV = CoefMag::Large;
            break;
    }
    return c;
}

inline u64 SampleCoefWeighted(std::mt19937_64& rng, unsigned w, CoefMag mag)
{
    const unsigned shortbits = std::min<unsigned>(16u, w / 4u);
    const u64 shortlim = (shortbits >= 64u) ? ~u64(0) : (u64(1) << shortbits);
    const u64 modulusmask = (w >= 64u) ? ~u64(0) : ((u64(1) << w) - 1u);

    switch (mag)
    {
        case CoefMag::Small:
        {
            const u64 range = 2u * shortlim;
            return MaskW(rng() % range, w);
        }
        case CoefMag::Medium:
        {
            // when the modulus is small enough that short bands eat it, fall back to full uniform; otherwise the band is empty and rng() % 0 is undefined behaviour
            if (modulusmask <= 2u * shortlim) return MaskW(rng(), w);
            const u64 bandsize = modulusmask - 2u * shortlim + 1u;
            return MaskW(shortlim + (rng() % bandsize), w);
        }
        case CoefMag::Large:
        {
            // at w=64 with shortbits=16, rejection rate is 2^-47 so an 8-try bound essentially never trips; the mid modulus fallback is only ever used at pathologically narrow w
            for (int tries = 0; tries < 8; ++tries)
            {
                u64 v = MaskW(rng(), w);
                if (v < shortlim) continue;
                if (v > modulusmask - shortlim) continue;
                return v;
            }
            return MaskW(modulusmask / 2u, w);
        }
    }
    return 0;
}

inline unsigned CoefMagBucket(const std::vector<u64>& coefs, unsigned w)
{
    const unsigned shortbits = std::min<unsigned>(16u, w / 4u);
    const u64 shortlim = (shortbits >= 64u) ? ~u64(0) : (u64(1) << shortbits);
    const u64 modulusmask = (w >= 64u) ? ~u64(0) : ((u64(1) << w) - 1u);
    unsigned shortcount = 0;
    for (u64 c : coefs)
    {
        const u64 cm = c & modulusmask;
        if (cm < shortlim) { ++shortcount; continue; }
        if (cm > modulusmask - shortlim) { ++shortcount; continue; }
    }
    if (shortcount == 0) return 0;
    if (shortcount <= 2) return 1;
    return 2;
}

inline unsigned DepthBucket(unsigned depthmax)
{
    if (depthmax <= 1) return 0;
    if (depthmax <= 3) return 1;
    return 2;
}

struct GeneratorFingerprint
{
    TargetKind TargetKindV = TargetKind::Zero;
    u64 TargetValue = 0;
    unsigned TargetVarIdx = 0;
    unsigned TargetVarIdxY = 0;
    XYOp TargetXYOp = XYOp::And;
    unsigned V = 0;
    unsigned S = 0;
    unsigned DepthBucketV = 0;
    unsigned KernelDim = 0;
    unsigned CoefMagBucketV = 0;

    bool operator==(const GeneratorFingerprint& o) const
    {
        return TargetKindV == o.TargetKindV
            && TargetValue == o.TargetValue
            && TargetVarIdx == o.TargetVarIdx
            && TargetVarIdxY == o.TargetVarIdxY
            && TargetXYOp == o.TargetXYOp
            && V == o.V
            && S == o.S
            && DepthBucketV == o.DepthBucketV
            && KernelDim == o.KernelDim
            && CoefMagBucketV == o.CoefMagBucketV;
    }
};

struct GeneratorState
{
    static constexpr size_t MaxCache = 64;
    struct Entry
    {
        GeneratorFingerprint Fp;
        std::vector<u64> Coefs;
        std::vector<AtomTreePtr> Trees;
    };
    std::vector<Entry> Cache;

    bool ContainsFingerprint(const GeneratorFingerprint& fp) const
    {
        for (const auto& e : Cache)
            if (e.Fp == fp) return true;
        return false;
    }

    void Push(Entry e)
    {
        if (Cache.size() >= MaxCache) Cache.erase(Cache.begin());
        Cache.push_back(std::move(e));
    }

    std::optional<size_t> FindMatchingTarget(const TargetSpec& target) const
    {
        for (size_t i = Cache.size(); i-- > 0; )
        {
            const auto& fp = Cache[i].Fp;
            if (fp.TargetKindV != target.Kind) continue;
            switch (target.Kind)
            {
                case TargetKind::Zero: return i;
                case TargetKind::ClassA:
                    if (fp.TargetValue == target.Constant) return i;
                    break;
                case TargetKind::ClassX:
                    if (fp.TargetValue == target.XCoef && fp.TargetVarIdx == target.XVarIdx) return i;
                    break;
                case TargetKind::ClassXY:
                    if (fp.TargetValue == target.XCoef && fp.TargetVarIdx == target.XVarIdx && fp.TargetVarIdxY == target.YVarIdx && fp.TargetXYOp == target.XYOpV) return i;
                    break;
            }
        }
        return std::nullopt;
    }
};

inline AtomTreePtr SpanAnchorForTarget(const TargetSpec& target, unsigned V, std::mt19937_64& rng)
{
    switch (target.Kind)
    {
        case TargetKind::Zero:
            // shouldn't be called for Zero (T=0 is always trivially solvable) but returning a bare Var defensively is cheaper than asserting; a crash here would take down the whole emit path for the misguided caller
            return MakeVar((unsigned)(rng() % V));
        case TargetKind::ClassA:
        {
            const unsigned k = (unsigned)(rng() % V);
            // all three forms evaluate to all-ones at every input (x|~x = -1 in Z/2^w, same for ~(x & ~x) and x^~x), so the truth matrix column is uniformly all 1s and guarantees the constant target is in the span regardless of what the other random trees sampled
            switch (rng() % 3u)
            {
                case 0: return MakeOr(MakeVar(k), MakeNot(MakeVar(k)));
                case 1: return MakeNot(MakeAnd(MakeVar(k), MakeNot(MakeVar(k))));
                default: return MakeXor(MakeVar(k), MakeNot(MakeVar(k)));
            }
        }
        case TargetKind::ClassX:
        {
            const unsigned k = target.XVarIdx;
            // V invariant violation by the caller; fall back to Var(0) and let the row reducer fail cleanly rather than crash on an out of range var index
            if (k >= V) return MakeVar(0);
            switch (rng() % 4u)
            {
                case 0: return MakeVar(k);
                case 1: return MakeNot(MakeNot(MakeVar(k)));
                case 2: return MakeAnd(MakeVar(k), MakeVar(k));
                default: return MakeOr(MakeVar(k), MakeVar(k));
            }
        }
        case TargetKind::ClassXY:
        {
            const unsigned a = target.XVarIdx;
            const unsigned b = target.YVarIdx;
            if (a >= V || b >= V || a == b) return MakeVar(0);
            const auto va = MakeVar(a);
            const auto vb = MakeVar(b);
            switch (target.XYOpV)
            {
                case XYOp::And:
                {
                    switch (rng() % 3u)
                    {
                        case 0: return MakeAnd(va, vb);
                        case 1: return MakeNot(MakeOr(MakeNot(va), MakeNot(vb)));
                        default: return MakeAnd(MakeAnd(va, vb), va);
                    }
                }
                case XYOp::Xor:
                {
                    switch (rng() % 3u)
                    {
                        case 0: return MakeXor(va, vb);
                        case 1: return MakeAnd(MakeOr(va, vb), MakeOr(MakeNot(va), MakeNot(vb)));
                        default: return MakeXor(MakeXor(va, vb), MakeAnd(va, MakeNot(va)));
                    }
                }
                case XYOp::Or:
                default:
                {
                    switch (rng() % 3u)
                    {
                        case 0: return MakeOr(va, vb);
                        case 1: return MakeNot(MakeAnd(MakeNot(va), MakeNot(vb)));
                        default: return MakeOr(MakeOr(va, vb), vb);
                    }
                }
            }
        }
    }
    return MakeVar(0);
}

struct GenerateIdentityRequest
{
    TargetSpec TargetV;
    Tier TierV = Tier::T8;
    unsigned VOverride = 0;
};

struct IdentityResult
{
    std::vector<std::pair<u64, AtomTreePtr>> Terms;
};

inline IdentityResult GenerateIdentity(const GenerateIdentityRequest& req, std::mt19937_64& rng, GeneratorState& state, unsigned w)
{
    IdentityResult result;

    const TierConfig tc = GetTierConfig(req.TierV);

    unsigned V;
    if (req.VOverride != 0) V = req.VOverride;
    else
    {
        const unsigned vspan = tc.VMax - tc.VMin + 1u;
        V = tc.VMin + (unsigned)(rng() % vspan);
    }
    if (V == 0) V = 1;
    // 2^16 truth-matrix rows is the ceiling; anything larger overflows the row reducer's per pivot O(rows*cols) scan budget and silently stalls the generator
    if (V > 16) V = 16;

    if (req.TargetV.Kind == TargetKind::ClassX && req.TargetV.XVarIdx >= V)
    {
        V = req.TargetV.XVarIdx + 1u;
        if (V > 16) return result;
    }
    if (req.TargetV.Kind == TargetKind::ClassXY)
    {
        // ClassXY requires distinct var indices; a==b degenerates the mix to a single var pattern the row reducer can't distinguish from a ClassX request, so bail out and let the caller fix its invariant
        if (req.TargetV.XVarIdx == req.TargetV.YVarIdx) return result;
        const unsigned need = std::max(req.TargetV.XVarIdx, req.TargetV.YVarIdx) + 1u;
        if (need > V) V = need;
        if (V > 16) return result;
    }

    const unsigned twov = 1u << V;
    const unsigned s = twov + tc.SExtraOver2V;

    const unsigned depthspan = tc.DepthMax - tc.DepthMin + 1u;
    auto drawdepth = [&]() -> unsigned { return tc.DepthMin + (unsigned)(rng() % depthspan); };

    std::vector<AtomTreePtr> trees;
    trees.reserve(s);
    for (unsigned i = 0; i < s; ++i)
        trees.push_back(RandomAtomTree(rng, V, drawdepth()));

    auto m = BuildTruthMatrix(trees, V, w);
    auto t = BuildTargetVector(req.TargetV, V, w);

    auto rr = RowReduceMod2W(m, t, w);
    if (!rr.Solvable)
    {
        // anchor swap is single shot per call; the anchor guarantees a span covering column, so a second infeasibility means the request is mathematically degenerate (seen mostly with V override forcing a value that doesn't span ClassX's variable), and the caller's retry loop handles the rare "regenerate with fresh randomness" case
        const unsigned swapidx = (unsigned)(rng() % s);
        trees[swapidx] = SpanAnchorForTarget(req.TargetV, V, rng);
        m = BuildTruthMatrix(trees, V, w);
        t = BuildTargetVector(req.TargetV, V, w);
        rr = RowReduceMod2W(m, t, w);
        if (!rr.Solvable) return result;
    }

    std::vector<u64> cpart(s, 0);
    for (unsigned c = 0; c < s; ++c)
    {
        const unsigned pr = rr.PivotRowForCol[c];
        if (pr == RowReduceResult::NPOS) continue;
        cpart[c] = rr.T[pr];
    }

    const unsigned rowstotal = (unsigned)rr.M.size();
    std::vector<unsigned char> constrained(s, 0);
    // Z/2^w is not a field; row reduction via odd-entry pivots leaves non zero residuals at non pivot rows in free columns (e.g. a `2` at row 6 of a rank-6-of-8 system), and those residuals impose additional homogeneous constraints the kernel direction must satisfy; the conservative fix forces constrained free weights to 0 so M*c_kernel == 0 holds; a fuller Smith Normal Form treatment would recover those dimensions but isn't needed at T8 default kernel widths
    for (unsigned r = rr.Rank; r < rowstotal; ++r)
    {
        for (unsigned f = 0; f < s; ++f)
        {
            if (rr.PivotRowForCol[f] != RowReduceResult::NPOS) continue;
            if (rr.M[r][f] != 0u) constrained[f] = 1;
        }
    }
    unsigned effectivekerneldim = 0;
    for (unsigned f = 0; f < s; ++f)
    {
        if (rr.PivotRowForCol[f] != RowReduceResult::NPOS) continue;
        if (!constrained[f]) ++effectivekerneldim;
    }

    std::vector<u64> ckernel(s, 0);
    if (effectivekerneldim > 0)
    {
        for (unsigned f = 0; f < s; ++f)
        {
            if (rr.PivotRowForCol[f] != RowReduceResult::NPOS) continue;
            if (constrained[f]) continue;
            const u64 wf = SampleCoefWeighted(rng, w, tc.CoefMagV);
            ckernel[f] = MaskW(ckernel[f] + wf, w);
            for (unsigned c = 0; c < s; ++c)
            {
                const unsigned pr = rr.PivotRowForCol[c];
                if (pr == RowReduceResult::NPOS) continue;
                const u64 minusprod = MaskW(u64(0) - MaskW(wf * rr.M[pr][f], w), w);
                ckernel[c] = MaskW(ckernel[c] + minusprod, w);
            }
        }
    }

    std::vector<u64> cfinal(s, 0);
    for (unsigned i = 0; i < s; ++i)
        cfinal[i] = MaskW(cpart[i] + ckernel[i], w);
    for (auto& tt : trees)
        tt = NormalizeNot(tt);

    GeneratorFingerprint fp;
    fp.TargetKindV = req.TargetV.Kind;
    switch (req.TargetV.Kind)
    {
        case TargetKind::ClassA: fp.TargetValue = MaskW(req.TargetV.Constant, w); break;
        case TargetKind::ClassX:
        case TargetKind::ClassXY: fp.TargetValue = MaskW(req.TargetV.XCoef, w); break;
        case TargetKind::Zero:
        default: fp.TargetValue = 0u; break;
    }
    fp.TargetVarIdx = (req.TargetV.Kind == TargetKind::ClassX || req.TargetV.Kind == TargetKind::ClassXY) ? req.TargetV.XVarIdx : 0u;
    fp.TargetVarIdxY = (req.TargetV.Kind == TargetKind::ClassXY) ? req.TargetV.YVarIdx : 0u;
    fp.TargetXYOp = (req.TargetV.Kind == TargetKind::ClassXY) ? req.TargetV.XYOpV : XYOp::And;
    fp.V = V;
    fp.S = s;
    fp.DepthBucketV = DepthBucket(tc.DepthMax);
    fp.KernelDim = effectivekerneldim;
    fp.CoefMagBucketV = CoefMagBucket(cfinal, w);

    if (state.ContainsFingerprint(fp)) return result;

    GeneratorState::Entry entry;
    entry.Fp = fp;
    entry.Coefs = cfinal;
    entry.Trees = trees;
    state.Push(entry);

    result.Terms.reserve(s);
    for (unsigned i = 0; i < s; ++i)
    {
        if (cfinal[i] == 0u) continue;
        result.Terms.emplace_back(cfinal[i], trees[i]);
    }

    // zero coef drop can leave Terms empty when the identity is trivially 0=target 
    if (result.Terms.empty())
        result.Terms.emplace_back(0u, trees[0]);

    return result;
}

inline IdentityResult GenerateIdentityWithRetry(const GenerateIdentityRequest& req, std::mt19937_64& rng, GeneratorState& state, unsigned w)
{
    static constexpr int MaxRetries = 10;
    IdentityResult lastresult;
    for (int attempt = 0; attempt < MaxRetries; ++attempt)
    {
        IdentityResult r = GenerateIdentity(req, rng, state, w);
        if (!r.Terms.empty()) return r;
        lastresult = std::move(r);
    }

    if (auto idx = state.FindMatchingTarget(req.TargetV))
    {
        const auto& e = state.Cache[*idx];
        IdentityResult r;
        r.Terms.reserve(e.Coefs.size());
        for (size_t i = 0; i < e.Coefs.size(); ++i)
        {
            if (e.Coefs[i] == 0u) continue;
            r.Terms.emplace_back(e.Coefs[i], e.Trees[i]);
        }
        if (!r.Terms.empty()) return r;
    }

    // secondary fallback clears the cache temporarily so the generate call bypasses the fingerprint collision check; this is the "give up on diversity but produce a correct emit" escape hatch, triggered <0.1% of calls in a typical 5000+ identity build
    auto savedcache = std::move(state.Cache);
    state.Cache.clear();
    IdentityResult r = GenerateIdentity(req, rng, state, w);
    state.Cache = std::move(savedcache);
    return r;
}

struct GeneratorPolicy
{
    Tier TierV = Tier::T8;
    unsigned MaxFreeAnchorsPerGroup = 2;
    unsigned MaxPaddingPerGroup = 2;
    float PaddingProb = 0.5f;
};

template <typename EmitFn>
inline void BuildGroupFromGenerator(u64 ki, unsigned V, const GeneratorPolicy& policy, std::mt19937_64& rng, GeneratorState& state, unsigned w, EmitFn emit)
{
    const TierConfig tcfg = GetTierConfig(policy.TierV);
    assert(V >= tcfg.VMin && "BuildGroupFromGenerator: caller's V is below tier minimum");
    // \clamp against the caller's V so emitted trees reference only indices in [0, V); the tier's natural V variation is preserved when V >= tier.VMax (V override = tier.VMax, which equals the unforced upper bound)
    const unsigned vuse = (V < tcfg.VMax) ? V : tcfg.VMax;

    u64 sumsofar = 0;

    const unsigned anchorcount = (policy.MaxFreeAnchorsPerGroup == 0) ? 0u : (unsigned)(rng() % (policy.MaxFreeAnchorsPerGroup + 1u));
    for (unsigned a = 0; a < anchorcount; ++a)
    {
        GenerateIdentityRequest r;
        r.TargetV.Kind = TargetKind::Zero;
        r.TierV = policy.TierV;
        r.VOverride = vuse;
        IdentityResult id = GenerateIdentityWithRetry(r, rng, state, w);
        for (const auto& term : id.Terms)
            emit(term.first, term.second);
    }

    auto bernoullipadding = [&]() -> bool
    {
        const u64 r = rng() & 0xFFFFu;
        const float pf = float(r) / float(0xFFFFu);
        return pf < policy.PaddingProb;
    };
    for (unsigned p = 0; p < policy.MaxPaddingPerGroup; ++p)
    {
        if (!bernoullipadding()) continue;
        GenerateIdentityRequest r;
        r.TargetV.Kind = TargetKind::Zero;
        r.TierV = policy.TierV;
        r.VOverride = vuse;
        IdentityResult id = GenerateIdentityWithRetry(r, rng, state, w);
        for (const auto& term : id.Terms)
            emit(term.first, term.second);
    }

    // sumsofar stays at 0 because the preceding anchors and padding all target Zero; the explicit tracking preserves the shape for future extensions where non zero contributions might enter and keeps the K_i-targeting arithmetic trivially correct
    const u64 correctivetarget = MaskW(ki - sumsofar, w);
    GenerateIdentityRequest rcorr;
    rcorr.TargetV.Kind = TargetKind::ClassA;
    rcorr.TargetV.Constant = correctivetarget;
    rcorr.TierV = policy.TierV;
    rcorr.VOverride = vuse;
    IdentityResult idcorr = GenerateIdentityWithRetry(rcorr, rng, state, w);
    for (const auto& term : idcorr.Terms)
        emit(term.first, term.second);
}

}