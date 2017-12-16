#include <stdexcept>

#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <thread>
#include <tuple>
#include <vector>

#include <pbbam/FastqReader.h>
#include <pbbam/Position.h>
#include <pbbam/QualityValues.h>

#include "Index.h"
#include "Mapping.h"

using namespace PacBio::minimap2;

using PacBio::BAM::FastqReader;
using PacBio::BAM::FastqSequence;
using PacBio::BAM::QualityValues;
using PacBio::BAM::Position;

#define mm_seq4_get(s, i)    ((s)[(i)>>3] >> (((i)&7)<<2) & 0xf)

void ReportCigar(const mm_reg1_t* const r, const int qlen, const int opt_flag, const Index& idx, const FastqSequence& record)
{
    constexpr const auto bases = "ACGT";
    constexpr const auto True = "True";
    constexpr const auto False = "False";

    uint32_t k, clip_len[2];
    clip_len[0] = r->rev? qlen - r->qe : r->qs;
    clip_len[1] = r->rev? r->qs : qlen - r->qe;
    const char clip_char = !(opt_flag & MM_F_SOFTCLIP) ? 'H' : 'S';  /* (sam_flag & 0x800) && */

    const std::string seq = record.Bases();
    const QualityValues quals = record.Qualities();
    const std::string refName{idx.idx_->seq[r->rid].name};
    const std::string prefix = "NA,NA," + refName;
    const Position refLen = idx.idx_->seq[r->rid].len;

    const auto refDatum = [&](const Position refPos) {
        return mm_seq4_get(idx.idx_->S, idx.idx_->seq[r->rid].offset + refPos);
    };

    const auto refBases = [&](const Position refPos, const int len) {
        std::ostringstream s;
        for (int i = 0; i < len; ++i)
            s << bases[refDatum(refPos + i)];
        return s.str();
    };

    const auto meanQV = [&](const Position qryPos, const int len) {
        double meanErr = 0.0;
        for (int i = 0; i < len; ++i)
            meanErr += std::pow(10.0, static_cast<double>(quals[qryPos + i]) / -10.0);
        meanErr /= len;
        return static_cast<int>(std::round(-10.0 * std::log10(meanErr)));
    };

    const char* inHP;
    int hpLen;
    char hpChr;
    const auto hpData = [&](const Position refPos, const std::string qryBases) {
        const auto d = refDatum(refPos);
        Position x, y;
        for (x = refPos; x > 0 && refDatum(x - 1) == d; --x);
        for (y = refPos + 1; y < refLen && refDatum(y) == d; ++y);
        const int l = y - x;
        const char* hp = False;
        if (qryBases.empty()) { if (l > 1) hp = True; }
        else if (qryBases.back() == bases[d]) hp = True;
        return std::make_tuple(hp, l, bases[d]);
    };

    Position qryPos = 0;
    Position refPos = r->rs;

    if (clip_len[0]) {
        const int len = clip_len[0];
        std::tie(inHP, hpLen, hpChr) = hpData(refPos, seq.substr(qryPos, len));
        const std::string atEnd = (refPos == 0 || refPos == refLen) ? True : False;
        std::cout << prefix << "," << refPos << "," << len << ",INDEL," << quals[qryPos] << ",";
        std::cout << atEnd << ",NA,NA,Insertion," << inHP << "," << seq.substr(qryPos, len) << "," << hpLen << "," << hpChr << std::endl;
        qryPos += len;
    }

    // "Movie,ZMW,Ref,Pos,Length,Type,QV,AtEnd,RefBP,AltBP,IndelType,InHP,Bases,HPLen,HPChar"
    for (k = 0; k < r->p->n_cigar; ++k) {
        const int len = r->p->cigar[k] >> 4;
        const std::string atEnd = (refPos == 0 || refPos == refLen) ? True : False;
        // MIDN=X
        switch (r->p->cigar[k] & 0xf) {
            case 0:  // MATCH
                qryPos += len;
                refPos += len;
                break;
            case 1:  // INSERTION
                std::tie(inHP, hpLen, hpChr) = hpData(refPos, seq.substr(qryPos, len));
                std::cout << prefix << "," << refPos << "," << len << ",INDEL," << meanQV(qryPos, len) << ",";
                std::cout << atEnd << ",NA,NA,Insertion," << inHP << "," << seq.substr(qryPos, len) << "," << hpLen << "," << hpChr << std::endl;
                qryPos += len;
                break;
            case 2:  // DELETION
                std::tie(inHP, hpLen, hpChr) = hpData(refPos, "");
                std::cout << prefix << "," << refPos << "," << len << ",INDEL," << meanQV(qryPos, 1) << ",";
                std::cout << atEnd << ",NA,NA,Deletion," << inHP << "," << refBases(refPos, len) << "," << hpLen << "," << hpChr << std::endl;
                refPos += len;
                break;
            case 3:  // SKIP
                refPos += len;
                std::cerr << "WARNING: encountered skip" << std::endl;
                break;
            case 4:  // MATCH
                qryPos += len;
                refPos += len;
                break;
            case 5:  // SUBSTITUTION
                std::cout << prefix << "," << refPos << "," << len << ",SNP," << meanQV(qryPos, len) << ",";
                std::cout << atEnd << "," << refBases(refPos, len) << "," << seq.substr(qryPos, len) << ",NA,NA,NA,NA,NA" << std::endl;
                qryPos += len;
                refPos += len;
                break;
            default:
                std::cerr << "ERROR: invalid cigar operation detected" << std::endl;
                exit(-1);
        }
    }

    if (clip_len[1]) {
        const int len = clip_len[0];
        std::tie(inHP, hpLen, hpChr) = hpData(refPos, seq.substr(qryPos, len));
        const std::string atEnd = (refPos == 0 || refPos == refLen) ? True : False;
        std::cout << prefix << "," << refPos << "," << len << ",INDEL," << quals[qryPos] << ",";
        std::cout << atEnd << ",NA,NA,Insertion," << inHP << "," << seq.substr(qryPos, len) << "," << hpLen << "," << hpChr << std::endl;
    }
}

int main(int argc, char* argv[])
{
    if (argc != 3) {
        std::cerr << "usage: rscheck REF.fasta CSS.fastq" << std::endl;
        return -1;
    }

    const std::string reffa(argv[1]), cssfq(argv[2]);

    IndexOptions idxOpts;
    // asm5 := io->flag = 0, io->k = 19, io->w = 19;
    idxOpts.opts_.flag = 0;
    idxOpts.opts_.k = 19;
    idxOpts.opts_.w = 19;

    MapOptions mapOpts;
    // asm5 := mo->a = 1; mo->b = 9, mo->q = 16, mo->q2 = 41, mo->e = 2, mo->e2 = 1, mo->zdrop = 200, mo->min_dp_max = 200, mo->best_n = 50;
    mapOpts.opts_.a = 1;
    mapOpts.opts_.b = 9;
    mapOpts.opts_.q = 16;
    mapOpts.opts_.q2 = 41;
    mapOpts.opts_.e = 2;
    mapOpts.opts_.e2 = 1;
    mapOpts.opts_.zdrop = 200;
    mapOpts.opts_.min_dp_max = 200;
    mapOpts.opts_.best_n = 50;

    {
        static constexpr const uint32_t nCigarMax = 65535;

        Index idx(reffa, idxOpts);
        mapOpts.Update(idx);

        ThreadBuffer tbuf;

        std::cout << "Movie,ZMW,Ref,Pos,Length,Type,QV,AtEnd,RefBP,AltBP,IndelType,InHP,Bases,HPLen,HPChar" << std::endl;

        FastqReader reader{cssfq};
        FastqSequence record;
        while (reader.GetNext(record)) {
            int numAlns;
            const auto seq = record.Bases();
            const int qlen = seq.length();
            auto alns = mm_map(idx.idx_, qlen, seq.c_str(), &numAlns, tbuf.tbuf_, &mapOpts.opts_, nullptr);
            for (int i = 0; i < numAlns; ++i) {
                auto aln = alns[i];
                // if no alignment, continue
                if (aln.p == nullptr) continue;
                if (aln.p->n_cigar > nCigarMax) {
                    std::cerr << "WARNING: hit cigar limit" << std::endl;
                    continue;
                }
                ReportCigar(&aln, qlen, mapOpts.opts_.flag, idx, record);
                /*
                const int32_t refId = aln.rid;
                const Position refStart = aln.rs;
                const Strand strand = aln.rev ? Strand::REVERSE : Strand::FORWARD;
                const Cigar cigar = RenderCigar(&aln, qlen, mapOpts.opts_.flag);
                const uint8_t mapq = aln.mapq;
                */
            }
            // cleanup
            for (int i = 0; i < numAlns; ++i) if (alns[i].p) free(alns[i].p);
            free(alns);
        }
    }

    return 0;
}