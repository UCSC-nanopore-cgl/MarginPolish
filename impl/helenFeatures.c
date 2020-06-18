//
// Created by tpesout on 3/29/19.
//

#ifdef _HDF5

#include "margin.h"
#include "htsIntegration.h"
#include "helenFeatures.h"
#include "ssw.h"
#include <hdf5.h>

#define TRUTH_ALN_LOG_LEVEL debug
#define TRUTH_ALN_IDENTITY_THRESHOLD .99
#define TRUTH_ALN_MIN_MATCHES 700

PoaFeatureSimpleWeight *PoaFeature_SimpleWeight_construct(int64_t refPos, int64_t insPos) {
    PoaFeatureSimpleWeight *feature = st_calloc(1, sizeof(PoaFeatureSimpleWeight));
    feature->refPosition = refPos;
    feature->insertPosition = insPos;
    feature->label = '\0';
    feature->nextInsert = NULL;
    return feature;
}

void PoaFeature_SimpleWeight_destruct(PoaFeatureSimpleWeight *feature) {
    if (feature->nextInsert != NULL) {
        PoaFeature_SimpleWeight_destruct(feature->nextInsert);
    }
    free(feature);
}

PoaFeatureSplitRleWeight *PoaFeature_SplitRleWeight_construct(int64_t refPos, int64_t insPos, int64_t rlPos,
                                                              int64_t maxRunLength) {
    PoaFeatureSplitRleWeight *feature = st_calloc(1, sizeof(PoaFeatureSplitRleWeight));
    feature->refPosition = refPos;
    feature->insertPosition = insPos;
    feature->runLengthPosition = rlPos;
    feature->labelChar = '\0';
    feature->labelRunLength = 0;
    feature->nextRunLength = NULL;
    feature->nextInsert = NULL;
    feature->maxRunLength = maxRunLength;
    feature->weights = st_calloc(((SYMBOL_NUMBER - 1) * (1 + maxRunLength) + 1) * 2, sizeof(double));
    return feature;
}

void PoaFeature_SplitRleWeight_destruct(PoaFeatureSplitRleWeight *feature) {
    if (feature->nextRunLength != NULL) {
        PoaFeature_SplitRleWeight_destruct(feature->nextRunLength);
    }
    if (feature->nextInsert != NULL) {
        PoaFeature_SplitRleWeight_destruct(feature->nextInsert);
    }
    free(feature->weights);
    free(feature);
}


PoaFeatureChannelRleWeight *PoaFeature_ChannelRleWeight_construct(int64_t refPos, int64_t insPos, int64_t rlPos,
                                                                  int64_t maxRunLength) {
    PoaFeatureChannelRleWeight *feature = st_calloc(1, sizeof(PoaFeatureChannelRleWeight));
    feature->refPosition = refPos;
    feature->insertPosition = insPos;
    feature->runLengthPosition = rlPos;
    feature->labelChar = '\0';
    feature->labelRunLength = 0;
    feature->nextRunLength = NULL;
    feature->nextInsert = NULL;
    feature->maxRunLength = maxRunLength;
    feature->nucleotideWeights = st_calloc((SYMBOL_NUMBER) * 2, sizeof(double));
    feature->runLengthWeights = st_calloc((SYMBOL_NUMBER - 1) * (1 + maxRunLength) * 2, sizeof(double));
    return feature;
}

void PoaFeature_ChannelRleWeight_destruct(PoaFeatureChannelRleWeight *feature) {
    if (feature->nextRunLength != NULL) {
        PoaFeature_ChannelRleWeight_destruct(feature->nextRunLength);
    }
    if (feature->nextInsert != NULL) {
        PoaFeature_ChannelRleWeight_destruct(feature->nextInsert);
    }
    free(feature->nucleotideWeights);
    free(feature->runLengthWeights);
    free(feature);
}

PoaFeatureDiploidRleWeight *PoaFeature_DiploidRleWeight_construct(int64_t refPos, int64_t insPos, int64_t rlPos,
        int64_t maxRunLength) {
    PoaFeatureDiploidRleWeight *feature = st_calloc(1, sizeof(PoaFeatureDiploidRleWeight));
    feature->refPosition = refPos;
    feature->insertPosition = insPos;
    feature->runLengthPosition = rlPos;
    feature->labelChar = '\0';
    feature->labelRunLength = 0;
    feature->nextRunLength = NULL;
    feature->nextInsert = NULL;
    feature->maxRunLength = maxRunLength;
    feature->weightsHOn = st_calloc(((SYMBOL_NUMBER - 1) * (1 + maxRunLength) + 1) * 2, sizeof(double));
    feature->weightsHOff = st_calloc(((SYMBOL_NUMBER - 1) * (1 + maxRunLength) + 1) * 2, sizeof(double));
    return feature;
}
void PoaFeature_DiploidRleWeight_destruct(PoaFeatureDiploidRleWeight *feature) {
    if (feature->nextRunLength != NULL) {
        PoaFeature_DiploidRleWeight_destruct(feature->nextRunLength);
    }
    if (feature->nextInsert != NULL) {
        PoaFeature_DiploidRleWeight_destruct(feature->nextInsert);
    }
    free(feature->weightsHOn);
    free(feature->weightsHOff);
    free(feature);
}


int PoaFeature_SimpleWeight_charIndex(Symbol character, bool forward) {
    int pos = character * 2 + (forward ? POS_STRAND_IDX : NEG_STRAND_IDX);
    assert(pos < POAFEATURE_SIMPLE_WEIGHT_TOTAL_SIZE);
    return pos;
}

int PoaFeature_SimpleWeight_gapIndex(bool forward) {
    int pos = POAFEATURE_SYMBOL_GAP_POS * 2 + (forward ? POS_STRAND_IDX : NEG_STRAND_IDX);
    assert(pos < POAFEATURE_SIMPLE_WEIGHT_TOTAL_SIZE);
    return pos;
}

int PoaFeature_SplitRleWeight_charIndex(int64_t maxRunLength, Symbol character, int64_t runLength, bool forward) {
    assert(runLength >= 0);
    assert(runLength <= maxRunLength);
    int pos = (character * ((int) maxRunLength + 1) + runLength) * 2 + (forward ? POS_STRAND_IDX : NEG_STRAND_IDX);
    return pos;
}

int PoaFeature_SplitRleWeight_gapIndex(int64_t maxRunLength, bool forward) {
    int pos = ((SYMBOL_NUMBER - 1) * ((int) maxRunLength + 1)) * 2 + (forward ? POS_STRAND_IDX : NEG_STRAND_IDX);
    return pos;
}

int PoaFeature_ChannelRleWeight_charNuclIndex(Symbol character, bool forward) {
    int pos = character * 2 + (forward ? POS_STRAND_IDX : NEG_STRAND_IDX);
    return pos;
}

int PoaFeature_ChannelRleWeight_gapNuclIndex(bool forward) {
    int pos = (SYMBOL_NUMBER_NO_N) * 2 + (forward ? POS_STRAND_IDX : NEG_STRAND_IDX);
    return pos;
}

int PoaFeature_ChannelRleWeight_charRLIndex(int64_t maxRunLength, Symbol character, int64_t runLength, bool forward) {
    assert(runLength >= 0);
    assert(runLength <= maxRunLength);
    int pos = (character * ((int) maxRunLength + 1) + runLength) * 2 + (forward ? POS_STRAND_IDX : NEG_STRAND_IDX);
    return pos;
}

int PoaFeature_DiploidRleWeight_charIndex(int64_t maxRunLength, Symbol character, int64_t runLength, bool forward) {
    return PoaFeature_SplitRleWeight_charIndex(maxRunLength, character, runLength, forward);
}
int PoaFeature_DiploidRleWeight_gapIndex(int64_t maxRunLength, bool forward) {
    return PoaFeature_SplitRleWeight_gapIndex(maxRunLength, forward);
}


void PoaFeature_handleHelenFeatures(
        // global params
        HelenFeatureType helenFeatureType, int64_t splitWeightMaxRunLength, void **helenHDF5Files,
        bool fullFeatureOutput, char *trueReferenceBam, RleString *originalReference, Params *params,

        // chunk params
        char *logIdentifier, int64_t chunkIdx, BamChunk *bamChunk, Poa *poa, stList *bamChunkReads,
        char *polishedConsensusString, RleString *polishedRleConsensus) {

    st_logInfo(">%s Performing feature generation for chunk.\n", logIdentifier);

    // get filename
    char *helenFeatureOutfileBase = NULL;
    switch (helenFeatureType) {
        case HFEAT_SIMPLE_WEIGHT:
            helenFeatureOutfileBase = stString_print("simpleWeight.C%05"PRId64".%s-%"PRId64"-%"PRId64,
                                                     chunkIdx, bamChunk->refSeqName,
                                                     bamChunk->chunkBoundaryStart, bamChunk->chunkBoundaryEnd);
            break;
        case HFEAT_SPLIT_RLE_WEIGHT:
            // name of folder, not of file
            helenFeatureOutfileBase = stString_print("splitRleWeight.C%05"PRId64".%s-%"PRId64"-%"PRId64,
                                                     chunkIdx, bamChunk->refSeqName,
                                                     bamChunk->chunkBoundaryStart, bamChunk->chunkBoundaryEnd);
            break;
        case HFEAT_CHANNEL_RLE_WEIGHT:
            // name of folder, not of file
            helenFeatureOutfileBase = stString_print("channelRleWeight.C%05"PRId64".%s-%"PRId64"-%"PRId64,
                                                     chunkIdx, bamChunk->refSeqName,
                                                     bamChunk->chunkBoundaryStart, bamChunk->chunkBoundaryEnd);
            break;
        default:
            st_errAbort("Unhandled HELEN feature type!\n");
    }

    // necessary to annotate poa with truth (if true reference BAM has been specified)
    stList *trueRefAlignment = NULL;
    RleString *trueRefRleString = NULL;
    bool validReferenceAlignment = FALSE;

    // get reference chunk
    if (trueReferenceBam != NULL) {
        // get alignment of true ref to assembly
        stList *trueRefReads = stList_construct3(0, (void (*)(void *)) bamChunkRead_destruct);
        stList *trueRefAligns = stList_construct3(0, (void (*)(void *)) stList_destruct);
        // construct new chunk
        BamChunk *trueRefBamChunk = bamChunk_copyConstruct(bamChunk);
        BamChunker *trueReferenceBamChunker = bamChunker_copyConstruct(bamChunk->parent);
        free(trueReferenceBamChunker->bamFile);
        trueReferenceBamChunker->bamFile = stString_copy(trueReferenceBam);
        trueRefBamChunk->parent = trueReferenceBamChunker;
        // get true ref as "read"
        uint32_t trueAlignmentCount = convertToReadsAndAlignments(trueRefBamChunk, originalReference, trueRefReads,
                                                                  trueRefAligns, NULL);

        if (trueAlignmentCount == 1) {
            BamChunkRead *trueRefRead = stList_get(trueRefReads, 0);
            stList *truthAlign = stList_get(trueRefAligns, 0);
            trueRefRleString = rleString_copy(trueRefRead->rleRead);

            // get consensus seq for original positions
            int64_t originalRefRleChunkStartPos = stIntTuple_get(stList_get(truthAlign, 0), 0);
            int64_t originalRefRleChunkEndPos = stIntTuple_get(stList_get(truthAlign, stList_length(truthAlign) - 1), 0);

            int64_t consensusAlnShift = -1;
            RleString *consensusRegion = getConsensusByEstimatedOriginalReferencePositions(originalReference,
                    polishedRleConsensus, trueRefRleString, originalRefRleChunkStartPos, originalRefRleChunkEndPos,
                    &consensusAlnShift);
            assert(consensusAlnShift);

            // get alignment
            double score_consensus, alignIdentity;
            trueRefAlignment = alignConsensusAndTruthRLEWithKmerAnchors(consensusRegion, trueRefRleString,
                    &score_consensus, params->polishParams);
            shiftAlignmentCoords(trueRefAlignment, 0, consensusAlnShift);
            rleString_destruct(consensusRegion);

            // quick fail
            if (stList_length(trueRefAlignment) <= TRUTH_ALN_MIN_MATCHES) {
                alignIdentity = -1;
            } else {
                // trim edges, calculate identity
                stList_removeInterval(trueRefAlignment, stList_length(trueRefAlignment) - 10, 10);
                stList_removeInterval(trueRefAlignment, 0, 10);
                alignIdentity = calculateAlignIdentity(polishedRleConsensus, trueRefRleString, trueRefAlignment);
            }

            // loggit
            if (st_getLogLevel() >= TRUTH_ALN_LOG_LEVEL) {
                char *consensusRaw = rleString_expand(polishedRleConsensus);
                char *truthRaw = rleString_expand(trueRefRleString);
                st_logInfo("\n");
                st_logInfo(" %s RAW Consensus (length %d):\n    %s\n", logIdentifier, strlen(consensusRaw), consensusRaw);
                st_logInfo(" %s RAW Truth (length %d):\n    %s\n", logIdentifier, strlen(truthRaw), truthRaw);
                st_logInfo(" %s Alignment of truth consensus:\n", logIdentifier);
                printMEAAlignment2(polishedRleConsensus, trueRefRleString, trueRefAlignment);
                st_logInfo("\n");
                free(consensusRaw);
                free(truthRaw);
            }

            if (alignIdentity < TRUTH_ALN_IDENTITY_THRESHOLD) {
                st_logInfo(" %s True reference alignment failed with %d matches and align identity %f\n", logIdentifier,
                           stList_length(trueRefAlignment), alignIdentity);
            } else {
                validReferenceAlignment = TRUE;
            }

        }

        stList_destruct(trueRefReads);
        stList_destruct(trueRefAligns);
        bamChunk_destruct(trueRefBamChunk);
        bamChunker_destruct(trueReferenceBamChunker);
    }

    // either write it, or note that we failed to find a valid reference alignment
    if (trueReferenceBam != NULL && !validReferenceAlignment) {
        st_logInfo(" %s No valid reference alignment was found, skipping HELEN feature output.\n", logIdentifier);
    } else {
        st_logInfo(" %s Writing HELEN features with filename base: %s\n", logIdentifier, helenFeatureOutfileBase);

        // write the actual features (type dependent)
        PoaFeature_writeHelenFeatures(helenFeatureType, poa, bamChunkReads, helenFeatureOutfileBase,
                                      bamChunk, trueRefAlignment, polishedRleConsensus, trueRefRleString,
                                      fullFeatureOutput, splitWeightMaxRunLength,
                                      (HelenFeatureHDF5FileInfo **) helenHDF5Files);

        // write the polished chunk in fasta format
        if (fullFeatureOutput) {
            char *chunkPolishedRefFilename = stString_print("%s.fa", helenFeatureOutfileBase);
            char *chunkPolishedRefContigName = stString_print("%s\t%"PRId64"\t%"PRId64"\t%s",
                                                              bamChunk->refSeqName,
                                                              bamChunk->chunkBoundaryStart,
                                                              bamChunk->chunkBoundaryEnd,
                                                              helenFeatureOutfileBase);
            FILE *chunkPolishedRefOutFh = fopen(chunkPolishedRefFilename, "w");
            fastaWrite(polishedConsensusString, chunkPolishedRefContigName, chunkPolishedRefOutFh);
            fclose(chunkPolishedRefOutFh);
            free(chunkPolishedRefFilename);
            free(chunkPolishedRefContigName);
        }
    }

    // cleanup
    free(helenFeatureOutfileBase);
    if (trueRefAlignment != NULL) stList_destruct(trueRefAlignment);
    if (trueRefRleString != NULL) rleString_destruct(trueRefRleString);
}

void getDiploidHaplotypeAlignmentsRAW2RLE(RleString *polishedRleConsensusH1, RleString *polishedRleConsensusH2,
                                          RleString *trueRefRleStringA, RleString *trueRefRleStringB,
                                          RleString **trueRefRleStringToHap1, RleString **trueRefRleStringToHap2,
                                          stList **trueRefAlignmentToHap1, stList **trueRefAlignmentToHap2,
                                          Params *params, char *logIdentifier) {

    char *polishedConsensusStringH1 = rleString_expand(polishedRleConsensusH1);
    char *polishedConsensusStringH2 = rleString_expand(polishedRleConsensusH2);
    char *trueRefExpandedA = rleString_expand(trueRefRleStringA);
    char *trueRefExpandedB = rleString_expand(trueRefRleStringB);

    // align all to all
    uint16_t score_trueA_polished1, score_trueA_polished2, score_trueB_polished1, score_trueB_polished2;

    stList *trueRefAlignmentRawSpace_Polished1TrueA = alignConsensusAndTruthSSW(
            polishedConsensusStringH1, trueRefExpandedA, &score_trueA_polished1);

    stList *trueRefAlignmentRawSpace_Polished1TrueB = alignConsensusAndTruthSSW(
            polishedConsensusStringH1, trueRefExpandedB, &score_trueB_polished1);

    stList *trueRefAlignmentRawSpace_Polished2TrueA= alignConsensusAndTruthSSW(
            polishedConsensusStringH2, trueRefExpandedA, &score_trueA_polished2);

    stList *trueRefAlignmentRawSpace_Polished2TrueB = alignConsensusAndTruthSSW(
            polishedConsensusStringH2, trueRefExpandedB, &score_trueB_polished2);

    // determine best alignment
    bool use_trueA_polished1 = score_trueA_polished1 + score_trueB_polished2 > score_trueA_polished2 + score_trueB_polished1;

    // convert to rleSpace if appropriate
    if (params->polishParams->useRunLengthEncoding) {

        // hap1
        *trueRefRleStringToHap1 = rleString_copy(use_trueA_polished1 ? trueRefRleStringA : trueRefRleStringB);

        uint64_t *polishedRleConsensus1_nonRleToRleCoordinateMap = rleString_getNonRleToRleCoordinateMap(
                polishedRleConsensusH1);
        uint64_t *trueRefRleStringToHap1_nonRleToRleCoordinateMap = rleString_getNonRleToRleCoordinateMap(
                *trueRefRleStringToHap1);

        *trueRefAlignmentToHap1 = runLengthEncodeAlignment(
                use_trueA_polished1 ? trueRefAlignmentRawSpace_Polished1TrueA : trueRefAlignmentRawSpace_Polished1TrueB,
                polishedRleConsensus1_nonRleToRleCoordinateMap, trueRefRleStringToHap1_nonRleToRleCoordinateMap);

        free(polishedRleConsensus1_nonRleToRleCoordinateMap);
        free(trueRefRleStringToHap1_nonRleToRleCoordinateMap);

        // hap2
        *trueRefRleStringToHap2 = rleString_copy(use_trueA_polished1 ? trueRefRleStringB : trueRefRleStringA);

        uint64_t *polishedRleConsensus2_nonRleToRleCoordinateMap = rleString_getNonRleToRleCoordinateMap(
                polishedRleConsensusH2);
        uint64_t *trueRefRleStringToHap2_nonRleToRleCoordinateMap = rleString_getNonRleToRleCoordinateMap(
                *trueRefRleStringToHap2);

        *trueRefAlignmentToHap2 = runLengthEncodeAlignment(
                use_trueA_polished1 ? trueRefAlignmentRawSpace_Polished2TrueB : trueRefAlignmentRawSpace_Polished2TrueA,
                polishedRleConsensus2_nonRleToRleCoordinateMap, trueRefRleStringToHap2_nonRleToRleCoordinateMap);

        free(polishedRleConsensus2_nonRleToRleCoordinateMap);
        free(trueRefRleStringToHap2_nonRleToRleCoordinateMap);

        //cleanup
        stList_destruct(trueRefAlignmentRawSpace_Polished1TrueA);
        stList_destruct(trueRefAlignmentRawSpace_Polished1TrueB);
        stList_destruct(trueRefAlignmentRawSpace_Polished2TrueA);
        stList_destruct(trueRefAlignmentRawSpace_Polished2TrueB);
    }

    // debugging
    if (st_getLogLevel() >= TRUTH_ALN_LOG_LEVEL) {
        st_logInfo(" %s Alignment of truth to Hap1:\n", logIdentifier);
        printMEAAlignment2(polishedRleConsensusH1, *trueRefRleStringToHap1, *trueRefAlignmentToHap1);
        char *consensusRawH1 = rleString_expand(polishedRleConsensusH1);
        char *truthRawToH1 = rleString_expand(*trueRefRleStringToHap1);
        st_logInfo(" %s RAW Consensus Seq H1:\n    %s\n", logIdentifier, consensusRawH1);
        st_logInfo(" %s RAW Truth Seq H1:\n    %s\n", logIdentifier, truthRawToH1);
        free(consensusRawH1);
        free(truthRawToH1);

        st_logInfo(" %s Alignment of truth to Hap2:\n", logIdentifier);
        printMEAAlignment2(polishedRleConsensusH2, *trueRefRleStringToHap2, *trueRefAlignmentToHap2);
        char *consensusRawH2 = rleString_expand(polishedRleConsensusH2);
        char *truthRawToH2 = rleString_expand(*trueRefRleStringToHap1);
        st_logInfo(" %s RAW Consensus Seq H2:\n    %s\n", logIdentifier, consensusRawH2);
        st_logInfo(" %s RAW Truth Seq H2:\n    %s\n", logIdentifier, truthRawToH2);
        free(consensusRawH2);
        free(truthRawToH2);
    }
}



void getDiploidHaplotypeAlignmentsRLE(RleString *polishedRleConsensusH1, RleString *polishedRleConsensusH2,
                                      RleString *trueRefRleStringA, RleString *trueRefRleStringB,
                                      RleString **trueRefRleStringToHap1, RleString **trueRefRleStringToHap2,
                                      stList **trueRefAlignmentToHap1, stList **trueRefAlignmentToHap2,
                                      Params *params, char *logIdentifier) {

    // align all to all
    double score_polished1_trueA, score_polished2_trueA, score_polished1_trueB, score_polished2_trueB;

    stList *trueRefAlignmentRLESpace_Polished1TrueA = alignConsensusAndTruthRLEWithKmerAnchors(
            polishedRleConsensusH1, trueRefRleStringA, &score_polished1_trueA, params->polishParams);

    stList *trueRefAlignmentRLESpace_Polished1TrueB = alignConsensusAndTruthRLEWithKmerAnchors(
            polishedRleConsensusH1, trueRefRleStringB, &score_polished1_trueB, params->polishParams);

    stList *trueRefAlignmentRLESpace_Polished2TrueA = alignConsensusAndTruthRLEWithKmerAnchors(
            polishedRleConsensusH2, trueRefRleStringA, &score_polished2_trueA, params->polishParams);

    stList *trueRefAlignmentRLESpace_Polished2TrueB = alignConsensusAndTruthRLEWithKmerAnchors(
            polishedRleConsensusH2, trueRefRleStringB, &score_polished2_trueB, params->polishParams);

    // cis or trans
    bool use_polished1_trueA = score_polished1_trueA + score_polished2_trueB > score_polished2_trueA + score_polished1_trueB;

    // hap1
    *trueRefRleStringToHap1 = rleString_copy(use_polished1_trueA ? trueRefRleStringA : trueRefRleStringB);
    *trueRefAlignmentToHap1 = use_polished1_trueA ?
                              trueRefAlignmentRLESpace_Polished1TrueA : trueRefAlignmentRLESpace_Polished1TrueB;

    // hap2
    *trueRefRleStringToHap2 = rleString_copy(use_polished1_trueA ? trueRefRleStringB : trueRefRleStringA);
    *trueRefAlignmentToHap2 = use_polished1_trueA ?
                              trueRefAlignmentRLESpace_Polished2TrueB : trueRefAlignmentRLESpace_Polished2TrueA;

    // debugging
    if (st_getLogLevel() >= TRUTH_ALN_LOG_LEVEL) {
        st_logInfo(" %s Alignment of truth to Hap1:\n", logIdentifier);
        printMEAAlignment2(polishedRleConsensusH1, *trueRefRleStringToHap1, *trueRefAlignmentToHap1);
        char *consensusRawH1 = rleString_expand(polishedRleConsensusH1);
        char *truthRawToH1 = rleString_expand(*trueRefRleStringToHap1);
        st_logInfo(" %s RAW Consensus Seq H1:\n    %s\n", logIdentifier, consensusRawH1);
        st_logInfo(" %s RAW Truth Seq H1:\n    %s\n", logIdentifier, truthRawToH1);
        free(consensusRawH1);
        free(truthRawToH1);

        st_logInfo(" %s Alignment of truth to Hap2:\n", logIdentifier);
        printMEAAlignment2(polishedRleConsensusH2, *trueRefRleStringToHap2, *trueRefAlignmentToHap2);
        char *consensusRawH2 = rleString_expand(polishedRleConsensusH2);
        char *truthRawToH2 = rleString_expand(*trueRefRleStringToHap1);
        st_logInfo(" %s RAW Consensus Seq H2:\n    %s\n", logIdentifier, consensusRawH2);
        st_logInfo(" %s RAW Truth Seq H2:\n    %s\n", logIdentifier, truthRawToH2);
        free(consensusRawH2);
        free(truthRawToH2);
    }

    //cleanup
    if (use_polished1_trueA) {
        stList_destruct(trueRefAlignmentRLESpace_Polished1TrueB);
        stList_destruct(trueRefAlignmentRLESpace_Polished2TrueA);
    } else {
        stList_destruct(trueRefAlignmentRLESpace_Polished1TrueA);
        stList_destruct(trueRefAlignmentRLESpace_Polished2TrueB);
    }
}

double calculateAlignIdentity(RleString *XRLE, RleString *YRLE, stList *alignedPairs) {
    if (stList_length(alignedPairs) == 0) {
        return 0.0;
    }
    char *X = XRLE->rleString;
    char *Y = YRLE->rleString;
    uint64_t *Xrl = XRLE->repeatCounts;
    uint64_t *Yrl = YRLE->repeatCounts;

    // stats to track
    int64_t matches = 0;
    int64_t mismatches = 0;
    int64_t xInserts = 0;
    int64_t yInserts = 0;

    // iterate over alignment
    stListIterator *alignmentItor = stList_getIterator(alignedPairs);
    stIntTuple *currAlign = stList_getNext(alignmentItor);
    int64_t posX = stIntTuple_get(currAlign, 0);
    int64_t posY = stIntTuple_get(currAlign, 1);

    while (TRUE) {
        if (currAlign == NULL) break;
        int64_t currAlignPosX = stIntTuple_get(currAlign, 0);
        int64_t currAlignPosY = stIntTuple_get(currAlign, 1);
        // Y gap / X insert
        if (posX < currAlignPosX) {
            posX++;
            xInserts += Xrl[posX];
        }

            // X gap / Y insert
        else if (posY < currAlignPosY) {
            posY++;
            yInserts += Yrl[posY];
        }

            // match
        else if (posX == currAlignPosX && posY == currAlignPosY) {
            if (tolower(X[posX]) == tolower(Y[posY])) {
                if (Xrl[posX] == Yrl[posY]) {
                    matches += Yrl[posY];
                } else if (Xrl[posX] > Yrl[posY]) {
                    matches += Yrl[posY];
                    mismatches += Xrl[posX] - Yrl[posY];
                } else {
                    matches += Xrl[posX];
                    mismatches += Yrl[posY] - Xrl[posX];
                }
            } else {
                if (Xrl[posX] == Yrl[posY]) {
                    mismatches += Yrl[posY];
                } else if (Xrl[posX] > Yrl[posY]) {
                    mismatches += Yrl[posY];
                    xInserts += Xrl[posX] - Yrl[posY];
                } else {
                    mismatches += Xrl[posX];
                    yInserts += Yrl[posY] - Xrl[posX];
                }
            }
            posX++;
            posY++;
            currAlign = stList_getNext(alignmentItor);
        }

            // should never happen
        else {
            assert(FALSE);
        }
    }

    stList_destructIterator(alignmentItor);

    return 1.0 * matches / (matches + mismatches + xInserts + yInserts);
}

bool alignToBestConsensus(RleString *trueRefRleString, RleString *polishedRleConsensusH1,
        RleString *polishedRleConsensusH2, int64_t consensusAlnShiftH1, int64_t consensusAlnShiftH2,
        stList *truthAlignmentsH1, stList *truthAlignmentsH2,
        stList *truthAlignmentDescriptors, Params *params, char *alignmentDesc, char *logIdentifier) {

    // for tracking success
    char *newAlignmentDesc;
    bool foundMatch = FALSE;

    // align to both haplotypes
    double score_consensusH1, score_consensusH2;
    stList *alignToH1 = alignConsensusAndTruthRLEWithKmerAnchors(polishedRleConsensusH1, trueRefRleString,
                                                                 &score_consensusH1, params->polishParams);
    stList *alignToH2 = alignConsensusAndTruthRLEWithKmerAnchors(polishedRleConsensusH2, trueRefRleString,
                                                                 &score_consensusH2, params->polishParams);

    // quick fail
    if (stList_length(alignToH1) <= TRUTH_ALN_MIN_MATCHES || stList_length(alignToH2) <= TRUTH_ALN_MIN_MATCHES) {
        // no good alignment
        newAlignmentDesc = stString_print("-0_%s", alignmentDesc);
        stList_destruct(alignToH1);
        stList_destruct(alignToH2);
        stList_append(truthAlignmentDescriptors, newAlignmentDesc);
        return FALSE;
    }

    // trim edges, calculate identity
    stList_removeInterval(alignToH1, stList_length(alignToH1) - 10, 10);
    stList_removeInterval(alignToH2, stList_length(alignToH2) - 10, 10);
    stList_removeInterval(alignToH1, 0, 10);
    stList_removeInterval(alignToH2, 0, 10);
    double alignIdentityH1 = calculateAlignIdentity(polishedRleConsensusH1, trueRefRleString, alignToH1);
    double alignIdentityH2 = calculateAlignIdentity(polishedRleConsensusH2, trueRefRleString, alignToH2);
    score_consensusH1 = stList_length(alignToH1) * alignIdentityH1;
    score_consensusH2 = stList_length(alignToH2) * alignIdentityH2;

    if (score_consensusH1 == score_consensusH2) {
        // no good alignment
        newAlignmentDesc = stString_print("-0_%s", alignmentDesc);
        st_logInfo(" %s Reference alignment for %s failed with identities H1:%f, H2:%f and scores H1:%f, H2:%f\n",
                logIdentifier, newAlignmentDesc, alignIdentityH1, alignIdentityH2, score_consensusH1, score_consensusH2);
        if (st_getLogLevel() >= TRUTH_ALN_LOG_LEVEL) {
            char *consensusRawH1 = rleString_expand(polishedRleConsensusH1);
            char *consensusRawH2 = rleString_expand(polishedRleConsensusH2);
            char *truthRaw = rleString_expand(trueRefRleString);
            st_logInfo("\n");
            st_logInfo(" %s RAW Consensus H1 (length %d):\n    %s\n", logIdentifier, strlen(consensusRawH1), consensusRawH1);
            st_logInfo(" %s RAW Consensus H2 (length %d):\n    %s\n", logIdentifier, strlen(consensusRawH2), consensusRawH2);
            st_logInfo(" %s RAW Truth %s (length %d):\n    %s\n", logIdentifier, newAlignmentDesc, strlen(truthRaw), truthRaw);
            st_logInfo(" %s Alignment of truth %s to Hap1:\n", logIdentifier, newAlignmentDesc);
            printMEAAlignment2(polishedRleConsensusH1, trueRefRleString, alignToH1);
            st_logInfo(" %s Alignment of truth %s to Hap2:\n", logIdentifier, newAlignmentDesc);
            printMEAAlignment2(polishedRleConsensusH2, trueRefRleString, alignToH2);
            st_logInfo("\n");
            free(consensusRawH2);
            free(consensusRawH1);
            free(truthRaw);
        }
        stList_destruct(alignToH1);
        stList_destruct(alignToH2);

    } else if (score_consensusH1 > score_consensusH2) {
        // alignment better to h1
        stList_destruct(alignToH2);

        if (score_consensusH1 < TRUTH_ALN_IDENTITY_THRESHOLD) {
            st_logInfo(" %s True reference alignment failed for -1_%s, align identity : %f\n",
                       logIdentifier, alignmentDesc, score_consensusH1);
            newAlignmentDesc = stString_print("-1_%s", alignmentDesc);
        } else {
            newAlignmentDesc = stString_print("+1_%s", alignmentDesc);
            int64_t startAlign = stIntTuple_get(stList_get(alignToH1, 0), 0);
            int64_t endAlign = stIntTuple_get(stList_get(alignToH1, stList_length(alignToH1) - 1), 0);
            stList_append(truthAlignmentsH1, HelenFeatureTruthAlignment_construct(startAlign, endAlign, alignToH1,
                                                                                  rleString_copy(trueRefRleString)));
            foundMatch = TRUE;
        }
        if (st_getLogLevel() >= TRUTH_ALN_LOG_LEVEL) {
            char *consensusRaw = rleString_expand(polishedRleConsensusH1);
            char *truthRaw = rleString_expand(trueRefRleString);
            st_logInfo("\n");
            st_logInfo(" %s RAW Consensus H1 (length %d):\n    %s\n", logIdentifier, strlen(consensusRaw), consensusRaw);
            st_logInfo(" %s RAW Truth %s (length %d):\n    %s\n", logIdentifier, newAlignmentDesc, strlen(truthRaw), truthRaw);
            st_logInfo(" %s Alignment of truth %s to Hap1:\n", logIdentifier, newAlignmentDesc);
            printMEAAlignment2(polishedRleConsensusH1, trueRefRleString, alignToH1);
            st_logInfo("\n");
            free(consensusRaw);
            free(truthRaw);
        }

        if (!foundMatch) {
            stList_destruct(alignToH1);
        } else {
            shiftAlignmentCoords(alignToH1, 0, consensusAlnShiftH1);
        }

    } else {
        // alignment better to h2
        stList_destruct(alignToH1);
        if (score_consensusH2 < TRUTH_ALN_IDENTITY_THRESHOLD) {
            st_logInfo(" %s True reference alignment failed for -2_%s, align identity : %f\n",
                       logIdentifier, alignmentDesc, score_consensusH2);

            newAlignmentDesc = stString_print("-2_%s", alignmentDesc);
        } else {
            newAlignmentDesc = stString_print("+2_%s", alignmentDesc);
            int64_t startAlign = stIntTuple_get(stList_get(alignToH2, 0), 0);
            int64_t endAlign = stIntTuple_get(stList_get(alignToH2, stList_length(alignToH2) - 1), 0);
            stList_append(truthAlignmentsH2, HelenFeatureTruthAlignment_construct(startAlign, endAlign, alignToH2,
                                                                                  rleString_copy(trueRefRleString)));
            foundMatch = TRUE;
        }

        if (st_getLogLevel() >= TRUTH_ALN_LOG_LEVEL) {
            char *consensusRaw = rleString_expand(polishedRleConsensusH2);
            char *truthRaw = rleString_expand(trueRefRleString);
            st_logInfo("\n");
            st_logInfo(" %s RAW Consensus H2 (length %d):\n    %s\n", logIdentifier, strlen(consensusRaw), consensusRaw);
            st_logInfo(" %s RAW Truth %s (length %d):\n    %s\n", logIdentifier, newAlignmentDesc, strlen(truthRaw), truthRaw);
            st_logInfo(" %s Alignment of truth %s to Hap2:\n", logIdentifier, newAlignmentDesc);
            printMEAAlignment2(polishedRleConsensusH2, trueRefRleString, alignToH2);
            st_logInfo("\n");
            free(truthRaw);
            free(consensusRaw);
        }

        if (!foundMatch) {
            stList_destruct(alignToH2);
        } else {
            shiftAlignmentCoords(alignToH2, 0, consensusAlnShiftH2);
        }
    }

    //TODO remove this hickey hackey business
    if (!foundMatch && score_consensusH1 < TRUTH_ALN_IDENTITY_THRESHOLD) {
        st_logInfo(" %s Attempting SSW truth alignment for %s.\n", logIdentifier, newAlignmentDesc);
        uint16_t scoreH1, scoreH2;

        char *rawConsensusH1 = rleString_expand(polishedRleConsensusH1);
        char *rawConsensusH2 = rleString_expand(polishedRleConsensusH2);
        char *rawTruth = rleString_expand(trueRefRleString);

        alignToH1 = alignConsensusAndTruthSSW(rawConsensusH1, rawTruth, &scoreH1);
        alignToH2 = alignConsensusAndTruthSSW(rawConsensusH2, rawTruth, &scoreH2);

        uint64_t *polishedRleConsensus_nonRleToRleCoordinateMapH1 = rleString_getNonRleToRleCoordinateMap(polishedRleConsensusH1);
        uint64_t *polishedRleConsensus_nonRleToRleCoordinateMapH2 = rleString_getNonRleToRleCoordinateMap(polishedRleConsensusH2);
        uint64_t *trueRefRleString_nonRleToRleCoordinateMap = rleString_getNonRleToRleCoordinateMap(trueRefRleString);
        stList *alignedPairsRawSSWToRLEH1 = runLengthEncodeAlignment(alignToH1,
                polishedRleConsensus_nonRleToRleCoordinateMapH1, trueRefRleString_nonRleToRleCoordinateMap);
        stList *alignedPairsRawSSWToRLEH2 = runLengthEncodeAlignment(alignToH2,
                polishedRleConsensus_nonRleToRleCoordinateMapH2, trueRefRleString_nonRleToRleCoordinateMap);
        stList_destruct(alignToH1);
        stList_destruct(alignToH2);
        alignToH1 = alignedPairsRawSSWToRLEH1;
        alignToH2 = alignedPairsRawSSWToRLEH2;

        free(rawConsensusH1);
        free(rawConsensusH2);
        free(rawTruth);
        free(polishedRleConsensus_nonRleToRleCoordinateMapH1);
        free(polishedRleConsensus_nonRleToRleCoordinateMapH2);
        free(trueRefRleString_nonRleToRleCoordinateMap);
        score_consensusH1 = calculateAlignIdentity(polishedRleConsensusH1, trueRefRleString, alignToH1);
        score_consensusH2 = calculateAlignIdentity(polishedRleConsensusH2, trueRefRleString, alignToH2);

        if (score_consensusH1 == score_consensusH2) {
            st_logInfo(" %s SSW reference alignment for %s failed with identities H1:%f, H2:%f\n", logIdentifier,
                       newAlignmentDesc, score_consensusH1, score_consensusH2);
            stList_destruct(alignToH1);
            stList_destruct(alignToH2);
        } else if (score_consensusH1 > score_consensusH2) {
            // alignment better to h1
            stList_destruct(alignToH2);

            if (score_consensusH2 >= TRUTH_ALN_IDENTITY_THRESHOLD) {
                free(newAlignmentDesc);
                newAlignmentDesc = stString_print("*1_%s", alignmentDesc);
                st_logInfo(" %s True reference alignment for %s succeeded with SSW after failing!\n", logIdentifier, newAlignmentDesc);
                int64_t startAlign = stIntTuple_get(stList_get(alignToH1, 0), 0);
                int64_t endAlign = stIntTuple_get(stList_get(alignToH1, stList_length(alignToH1) - 1), 0);
                stList_append(truthAlignmentsH1, HelenFeatureTruthAlignment_construct(startAlign, endAlign, alignToH1,
                                                                                      rleString_copy(trueRefRleString)));
                foundMatch = TRUE;
            }
            if (st_getLogLevel() >= TRUTH_ALN_LOG_LEVEL) {
                char *truthRaw = rleString_expand(trueRefRleString);
                st_logInfo("\n");
                st_logInfo(" %s SSW Alignment of truth %s to Hap1:\n", logIdentifier, newAlignmentDesc);
                printMEAAlignment2(polishedRleConsensusH1, trueRefRleString, alignToH1);
                st_logInfo("\n");
                free(truthRaw);
            }
            if (!foundMatch) {
                stList_destruct(alignToH1);
            } else {
                shiftAlignmentCoords(alignToH1, 0, consensusAlnShiftH1);
            }
        } else {
            // alignment better to h2
            stList_destruct(alignToH1);

            if (score_consensusH2 >= TRUTH_ALN_IDENTITY_THRESHOLD) {
                free(newAlignmentDesc);
                newAlignmentDesc = stString_print("*2_%s", alignmentDesc);
                st_logInfo(" %s True reference alignment for %s succeeded with SSW after failing!\n", logIdentifier, newAlignmentDesc);
                int64_t startAlign = stIntTuple_get(stList_get(alignToH2, 0), 0);
                int64_t endAlign = stIntTuple_get(stList_get(alignToH2, stList_length(alignToH2) - 1), 0);
                stList_append(truthAlignmentsH2, HelenFeatureTruthAlignment_construct(startAlign, endAlign, alignToH2,
                                                                                      rleString_copy(trueRefRleString)));
                foundMatch = TRUE;
            }
            if (st_getLogLevel() >= TRUTH_ALN_LOG_LEVEL) {
                char *truthRaw = rleString_expand(trueRefRleString);
                st_logInfo("\n");
                st_logInfo(" %s SSW Alignment of truth %s to Hap2:\n", logIdentifier, newAlignmentDesc);
                printMEAAlignment2(polishedRleConsensusH2, trueRefRleString, alignToH2);
                st_logInfo("\n");
                free(truthRaw);
            }
            if (!foundMatch) {
                stList_destruct(alignToH2);
            } else {
                shiftAlignmentCoords(alignToH2, 0, consensusAlnShiftH2);
            }
        }
    }

    // save to list
    stList_append(truthAlignmentDescriptors, newAlignmentDesc);
    return foundMatch;
}

RleString *getConsensusByEstimatedOriginalReferencePositions(RleString *originalReference, RleString *consensus,
        RleString *trueRefRleString, int64_t originalRefRleChunkStartPos, int64_t originalRefRleChunkEndPos,
        int64_t *rleEstimatedConsensusStartPos) {

    // get reference coords
    uint64_t *originalReferenceRLEMap = rleString_getRleToNonRleCoordinateMap(originalReference);
    int64_t originalRefRawChunkStartPos = originalReferenceRLEMap[originalRefRleChunkStartPos];
    int64_t originalRefRawChunkEndPos = originalReferenceRLEMap[originalRefRleChunkEndPos];
    int64_t originalRefRawEstStartPos = originalRefRawChunkStartPos * consensus->nonRleLength / originalReference->nonRleLength;
    int64_t originalRefRawEstEndPos = originalRefRawChunkEndPos * consensus->nonRleLength / originalReference->nonRleLength;

    // get estimated positions
    *rleEstimatedConsensusStartPos = -1;
    int64_t rleEstimatedConsensusEndPos = -1;
    int64_t pos = 0;
    for (int64_t i = 0; i < consensus->nonRleLength; i++) {
        if (pos <= originalRefRawEstStartPos) {
            *rleEstimatedConsensusStartPos = i;
        }
        if (pos >= originalRefRawEstEndPos) {
            rleEstimatedConsensusEndPos = i;
            break;
        }
        pos += consensus->repeatCounts[i];
    }
    if (rleEstimatedConsensusEndPos < 0) {
        rleEstimatedConsensusEndPos = consensus->length;
    }

    //cleanup
    free(originalReferenceRLEMap);

    // sanity check
    assert(*rleEstimatedConsensusStartPos >= 0);
    assert(rleEstimatedConsensusEndPos >= *rleEstimatedConsensusStartPos);

    // get consensus
    RleString *truncatedConsensus = rleString_copySubstring(consensus, (uint64_t) *rleEstimatedConsensusStartPos,
            (uint64_t) rleEstimatedConsensusEndPos - *rleEstimatedConsensusStartPos);

    // loggit
    if (st_getLogLevel() >= debug/*TRUTH_ALN_LOG_LEVEL*/) {
        char *logIdentifier = getLogIdentifier();
        st_logInfo(" %s Getting aligned region of consensus:\n", logIdentifier);
        st_logInfo("          T  rle_start:   %9"PRId64"   raw_chunk_start:%9"PRId64"\n",
                   originalRefRleChunkStartPos, originalRefRawChunkStartPos);
        st_logInfo("          T  rle_end:     %9"PRId64"   raw_chunk_end:  %9"PRId64"\n",
                   originalRefRleChunkEndPos, originalRefRawChunkEndPos);
        st_logInfo("          T  rle_est_len: %9"PRId64"   raw_est_len:    %9"PRId64"\n",
                   originalRefRleChunkEndPos - originalRefRleChunkStartPos,
                   originalRefRawChunkEndPos - originalRefRawChunkStartPos);
        st_logInfo("          T  rle_len:     %9"PRIu64"   raw_len:        %9"PRIu64"\n", trueRefRleString->length,
                   trueRefRleString->nonRleLength);
        st_logInfo("          C  rle_len:     %9"PRIu64"   raw_len:        %9"PRIu64"\n", truncatedConsensus->length,
                   truncatedConsensus->nonRleLength);
        free(logIdentifier);
    }

    return truncatedConsensus;
}

void shiftAlignmentCoords(stList *alignedPairs, int64_t tupleIdx, int64_t shift) {
    int64_t *tuple = NULL;
    for (int64_t i = 0; i < stList_length(alignedPairs); i++) {
        tuple = stList_get(alignedPairs, i);
        tuple[tupleIdx + 1] += shift;
    }
}

bool PoaFeature_handleDiploidHelenTruthAlignment(char *trueReferenceBamA, char *trueReferenceBamB, BamChunk *bamChunk,
        RleString *originalReference, Poa *poaH1, Poa *poaH2,
        RleString *polishedRleConsensusH1, RleString *polishedRleConsensusH2,
        stList *finalAlignmentsToH1, stList *finalAlignmentsToH2,
        Params *params, char *logIdentifier) {

    // get alignment of true ref to assembly
    stList *truthReadsA = stList_construct3(0, (void (*)(void *)) bamChunkRead_destruct);
    stList *truthReadsB = stList_construct3(0, (void (*)(void *)) bamChunkRead_destruct);
    stList *truthOriginalRefAlignA = stList_construct3(0, (void (*)(void *)) stList_destruct);
    stList *truthOriginalRefAlignB = stList_construct3(0, (void (*)(void *)) stList_destruct);
    // construct new chunk
    BamChunk *trueRefBamChunk = bamChunk_copyConstruct(bamChunk);
    BamChunker *trueReferenceBamChunker = bamChunker_copyConstruct(bamChunk->parent);
    free(trueReferenceBamChunker->bamFile);
    trueRefBamChunk->parent = trueReferenceBamChunker;
    // get true ref as reads
    trueReferenceBamChunker->bamFile = trueReferenceBamA;
    convertToReadsAndAlignments(trueRefBamChunk, originalReference, truthReadsA, truthOriginalRefAlignA, NULL);
    trueReferenceBamChunker->bamFile = trueReferenceBamB;
    convertToReadsAndAlignments(trueRefBamChunk, originalReference, truthReadsB, truthOriginalRefAlignB, NULL);
    trueReferenceBamChunker->bamFile = NULL;
    // tracking
    stList *truthAlignmentDescriptors = stList_construct3(0, free);
    stList *truthAlignmentsH1 = stList_construct();
    stList *truthAlignmentsH2 = stList_construct();

    //logging
    if (st_getLogLevel() >= TRUTH_ALN_LOG_LEVEL) {
        char *consensusRawH1 = rleString_expand(polishedRleConsensusH1);
        char *consensusRawH2 = rleString_expand(polishedRleConsensusH2);
        st_logInfo(" %s RAW Consensus Seq H1:\n    %s\n", logIdentifier, consensusRawH1);
        st_logInfo(" %s RAW Consensus Seq H2:\n    %s\n", logIdentifier, consensusRawH1);
        free(consensusRawH1);
        free(consensusRawH2);
    }

    // align truth reads A
    for (int64_t i = 0; i < stList_length(truthReadsA); i++) {
        BamChunkRead *trueRefBamChunkReadA = stList_get(truthReadsA, i);
        RleString *trueRefRleStringA = trueRefBamChunkReadA->rleRead;
        stList *truthAlign = stList_get(truthOriginalRefAlignA, i);

        // for logging
        char *alignmentDesc = stList_length(truthAlign) == 0 ? stString_copy("A,0-0,0") :
                              stString_print("A,%"PRId64"-%"PRId64",%"PRId64, stIntTuple_get(stList_get(truthAlign, 0), 0),
                                             stIntTuple_get(stList_get(truthAlign, stList_length(truthAlign) - 1), 0),
                                             trueRefRleStringA->length);

        // get consensus seq for original positions
        int64_t originalRefRleChunkStartPos = stIntTuple_get(stList_get(truthAlign, 0), 0);
        int64_t originalRefRleChunkEndPos = stIntTuple_get(stList_get(truthAlign, stList_length(truthAlign) - 1), 0);
        int64_t consensusH1StartPos = 0;
        int64_t consensusH2StartPos = 0;
        RleString *consensusRegionH1 = getConsensusByEstimatedOriginalReferencePositions(originalReference,
                polishedRleConsensusH1, trueRefRleStringA, originalRefRleChunkStartPos, originalRefRleChunkEndPos,
                &consensusH1StartPos);
        RleString *consensusRegionH2 = getConsensusByEstimatedOriginalReferencePositions(originalReference,
                polishedRleConsensusH2, trueRefRleStringA, originalRefRleChunkStartPos, originalRefRleChunkEndPos,
                &consensusH2StartPos);
        assert(consensusH1StartPos >= 0 && consensusH2StartPos >= 0);


        // save alignment
        alignToBestConsensus(trueRefRleStringA, consensusRegionH1, consensusRegionH2, consensusH1StartPos, consensusH2StartPos,
                             truthAlignmentsH1, truthAlignmentsH2, truthAlignmentDescriptors, params, alignmentDesc, logIdentifier);

        free(alignmentDesc); //copied in a2BC before saving
        rleString_destruct(consensusRegionH1);
        rleString_destruct(consensusRegionH2);
    }

    // align truthReads B
    for (int64_t i = 0; i < stList_length(truthReadsB); i++) {
        BamChunkRead *trueRefBamChunkReadB = stList_get(truthReadsB, i);
        RleString *trueRefRleStringB = trueRefBamChunkReadB->rleRead;
        stList *truthAlign = stList_get(truthOriginalRefAlignB, i);

        // for logging
        char *alignmentDesc = stList_length(truthAlign) == 0 ? stString_copy("B,0-0,0") :
                              stString_print("B,%"PRId64"-%"PRId64",%"PRId64, stIntTuple_get(stList_get(truthAlign, 0), 0),
                                             stIntTuple_get(stList_get(truthAlign, stList_length(truthAlign) - 1), 0),
                                             trueRefRleStringB->length);

        // get consensus seq for original positions
        int64_t originalRefRleChunkStartPos = stIntTuple_get(stList_get(truthAlign, 0), 0);
        int64_t originalRefRleChunkEndPos = stIntTuple_get(stList_get(truthAlign, stList_length(truthAlign) - 1), 0);
        int64_t consensusH1StartPos = 0;
        int64_t consensusH2StartPos = 0;
        RleString *consensusRegionH1 = getConsensusByEstimatedOriginalReferencePositions(originalReference,
                polishedRleConsensusH1, trueRefRleStringB, originalRefRleChunkStartPos, originalRefRleChunkEndPos,
                &consensusH1StartPos);
        RleString *consensusRegionH2 = getConsensusByEstimatedOriginalReferencePositions(originalReference,
                polishedRleConsensusH2, trueRefRleStringB, originalRefRleChunkStartPos, originalRefRleChunkEndPos,
                &consensusH2StartPos);
        assert(consensusH1StartPos >= 0 && consensusH2StartPos >= 0);

        // save alignment
        alignToBestConsensus(trueRefRleStringB, consensusRegionH1, consensusRegionH2, consensusH1StartPos, consensusH2StartPos,
                             truthAlignmentsH1, truthAlignmentsH2, truthAlignmentDescriptors, params, alignmentDesc, logIdentifier);

        free(alignmentDesc); //copied in a2BC before saving
        rleString_destruct(consensusRegionH1);
        rleString_destruct(consensusRegionH2);
    }

    // loggit
    char *allAlignDescriptors = stString_join2("; ", truthAlignmentDescriptors);
    st_logInfo(" %s Found total of %"PRId64" truth sequences: %s\n", logIdentifier, stList_length(truthAlignmentDescriptors),
            allAlignDescriptors);
    free(allAlignDescriptors);
    stList_destruct(truthAlignmentDescriptors);

    // did we find a viable alignment
    bool foundValidAlignment = FALSE;

    //filter for non-overlap
    stList_sort(truthAlignmentsH1, (int (*)(const void *, const void *))HelenFeatureTruthAlignment_cmp);
    for (int64_t i = 0; i < stList_length(truthAlignmentsH1); i++) {
        HelenFeatureTruthAlignment *hfta = stList_get(truthAlignmentsH1, i);
        if (stList_length(finalAlignmentsToH1) == 0 ||
            ((HelenFeatureTruthAlignment*)stList_peek(finalAlignmentsToH1))->endPosExcl <= hfta->startPosIncl) {
            stList_append(finalAlignmentsToH1, hfta);
            foundValidAlignment = TRUE;
        } else {
            HelenFeatureTruthAlignment_destruct(hfta);
        }
    }

    stList_sort(truthAlignmentsH2, (int (*)(const void *, const void *))HelenFeatureTruthAlignment_cmp);
    for (int64_t i = 0; i < stList_length(truthAlignmentsH2); i++) {
        HelenFeatureTruthAlignment *hfta = stList_get(truthAlignmentsH2, i);
        if (stList_length(finalAlignmentsToH2) == 0 ||
            ((HelenFeatureTruthAlignment*)stList_peek(finalAlignmentsToH2))->endPosExcl <= hfta->startPosIncl) {
            stList_append(finalAlignmentsToH2, hfta);
            foundValidAlignment = TRUE;
        } else {
            HelenFeatureTruthAlignment_destruct(hfta);
        }
    }

    //cleanump
    stList_destruct(truthAlignmentsH1);
    stList_destruct(truthAlignmentsH2);
    stList_destruct(truthReadsA);
    stList_destruct(truthOriginalRefAlignA);
    stList_destruct(truthReadsB);
    stList_destruct(truthOriginalRefAlignB);
    bamChunk_destruct(trueRefBamChunk);
    bamChunker_destruct(trueReferenceBamChunker);

    return foundValidAlignment;
}

void PoaFeature_handleDiploidHelenFeatures(
        // global params
        HelenFeatureType helenFeatureType,
        int64_t splitWeightMaxRunLength, void **helenHDF5Files, bool fullFeatureOutput,
        char *trueReferenceBamA, char *trueReferenceBamB, Params *params,

        // chunk params
        char *logIdentifier, int64_t chunkIdx, BamChunk *bamChunk, stList *bamChunkReads, Poa *poaH1, Poa *poaH2,
        stSet *readsInH1, stSet *readsInH2, RleString *polishedRleConsensusH1, RleString *polishedRleConsensusH2,
        RleString *originalReference) {

    st_logInfo(">%s Performing diploid feature generation for chunk.\n", logIdentifier);

    // get filename
    char *helenFeatureOutfileBase = NULL;
    switch (helenFeatureType) {
        case HFEAT_DIPLOID_RLE_WEIGHT:
            helenFeatureOutfileBase = stString_print("diploidRleWeight.C%05"PRId64".%s-%"PRId64"-%"PRId64,
                                                     chunkIdx, bamChunk->refSeqName,
                                                     bamChunk->chunkBoundaryStart, bamChunk->chunkBoundaryEnd);
            break;
        default:
            st_errAbort("Unhandled HELEN feature type!\n");
    }

    // annotate poa with truth (if true reference BAM has been specified)
    bool validReferenceAlignment = FALSE;
    stList *truthAlignmentsToHap1 = NULL;
    stList *truthAlignmentsToHap2 = NULL;

    // get reference chunk
    if (trueReferenceBamA != NULL && trueReferenceBamB != NULL) {
        truthAlignmentsToHap1 = stList_construct3(0, (void(*)(void*))HelenFeatureTruthAlignment_destruct);
        truthAlignmentsToHap2 = stList_construct3(0, (void(*)(void*))HelenFeatureTruthAlignment_destruct);

        // assign truth reads to hap1 and hap2: alignmentsToHapX's are populated now
        validReferenceAlignment = PoaFeature_handleDiploidHelenTruthAlignment(trueReferenceBamA, trueReferenceBamB,
                bamChunk, originalReference, poaH1, poaH2, polishedRleConsensusH1, polishedRleConsensusH2,
                truthAlignmentsToHap1, truthAlignmentsToHap2, params, logIdentifier);

    }

    // either write it, or note that we failed to find a valid reference alignment
    if (trueReferenceBamA != NULL && trueReferenceBamB != NULL && !validReferenceAlignment) {
        st_logInfo(" %s No valid reference alignment was found, skipping HELEN feature output.\n", logIdentifier);
    } else {
        st_logInfo(" %s Writing HELEN features with filename base: %s\n", logIdentifier, helenFeatureOutfileBase);

        // write the actual features (type dependent)
        stSet *readIdsInH1 = bamChunkRead_to_readName(readsInH1);
        stSet *readIdsInH2 = bamChunkRead_to_readName(readsInH2);
        PoaFeature_writeDiploidHelenFeatures(helenFeatureType, bamChunkReads, helenFeatureOutfileBase,
                                             bamChunk, poaH1, poaH2, readIdsInH1, readIdsInH2, truthAlignmentsToHap1,
                                             truthAlignmentsToHap2, splitWeightMaxRunLength,
                                             (HelenFeatureHDF5FileInfo **) helenHDF5Files);
        stSet_destruct(readIdsInH1);
        stSet_destruct(readIdsInH2);
    }

    // cleanup
    free(helenFeatureOutfileBase);
    if (truthAlignmentsToHap1 != NULL) stList_destruct(truthAlignmentsToHap1);
    if (truthAlignmentsToHap2 != NULL) stList_destruct(truthAlignmentsToHap2);
}

stList *PoaFeature_getSimpleWeightFeatures(Poa *poa, stList *bamChunkReads) {

    // initialize feature list
    stList *featureList = stList_construct3(0, (void (*)(void *)) PoaFeature_SimpleWeight_destruct);
    for (int64_t i = 1; i < stList_length(poa->nodes); i++) {
        stList_append(featureList, PoaFeature_SimpleWeight_construct(i - 1, 0));
    }

    // for logging (of errors)
    char *logIdentifier = getLogIdentifier();

    // iterate over all positions
    for (int64_t i = 0; i < stList_length(featureList); i++) {

        // get feature and node
        PoaFeatureSimpleWeight *feature = stList_get(featureList, i);
        PoaNode *node = stList_get(poa->nodes, i +
                                               1); //skip the first poa node, as it's always an 'N', so featureIdx and poaIdx are off by one

        // examine each observation
        stList *observations = node->observations;
        for (int64_t o = 0; o < stList_length(observations); o++) {

            // get data
            PoaBaseObservation *observation = stList_get(observations, o);
            BamChunkRead *bamChunkRead = stList_get(bamChunkReads, observation->readNo);
            Symbol symbol = poa->alphabet->convertCharToSymbol(bamChunkRead->rleRead->rleString[observation->offset]);
            bool forward = bamChunkRead->forwardStrand;

            // save weight
            feature->weights[PoaFeature_SimpleWeight_charIndex(symbol, forward)] += observation->weight;
        }

        // Deletes
        if (stList_length(node->deletes) > 0) {

            // iterate over all deletes
            for (int64_t j = 0; j < stList_length(node->deletes); j++) {
                PoaDelete *delete = stList_get(node->deletes, j);

                // Deletes start AFTER the current position, need to add counts/weights to nodes after the current node
                for (int64_t k = 1; k < delete->length; k++) {
                    if (i + k >= stList_length(poa->nodes)) {
                        st_logInfo(" %s Encountered DELETE that occurs after end of POA!\n", logIdentifier);
                        break;
                    }
                    PoaFeatureSimpleWeight *delFeature = stList_get(featureList, i + k);

                    delFeature->weights[PoaFeature_SimpleWeight_gapIndex(TRUE)] += delete->weightForwardStrand;
                    delFeature->weights[PoaFeature_SimpleWeight_gapIndex(FALSE)] += delete->weightReverseStrand;
                }
            }
        }

        // Inserts
        if (stList_length(node->inserts) > 0) {

            // iterate over all inserts
            for (int64_t j = 0; j < stList_length(node->inserts); j++) {
                PoaInsert *insert = stList_get(node->inserts, j);

                // get feature iterator
                PoaFeatureSimpleWeight *prevFeature = feature;
                for (int64_t k = 0; k < strlen(insert->insert->rleString); k++) {
                    // get current feature (or create if necessary)
                    PoaFeatureSimpleWeight *currFeature = prevFeature->nextInsert;
                    if (currFeature == NULL) {
                        currFeature = PoaFeature_SimpleWeight_construct(i, k + 1);
                        prevFeature->nextInsert = currFeature;
                    }

                    Symbol c = poa->alphabet->convertCharToSymbol(insert->insert->rleString[k]);

                    // add weights
                    currFeature->weights[PoaFeature_SimpleWeight_charIndex(c, TRUE)] += insert->weightForwardStrand;
                    currFeature->weights[PoaFeature_SimpleWeight_charIndex(c, FALSE)] += insert->weightReverseStrand;

                    // iterate
                    prevFeature = currFeature;
                }
            }
        }
    }

    free(logIdentifier);
    return featureList;
}


void poa_addSplitRunLengthFeaturesForObservations(Poa *poa, PoaFeatureSplitRleWeight *baseFeature, stList *observations,
                                                  stList *bamChunkReads, const int64_t maxRunLength,
                                                  int64_t observationOffset) {


    PoaFeatureSplitRleWeight *currFeature = baseFeature;
    int64_t currentRunLengthIndex = 0;
    bool beforeMaxObservedRunLength = TRUE;

    while (beforeMaxObservedRunLength) {
        beforeMaxObservedRunLength = FALSE;

        // examine each observation
        for (int64_t i = 0; i < stList_length(observations); i++) {

            // get data
            PoaBaseObservation *observation = stList_get(observations, i);
            BamChunkRead *bamChunkRead = stList_get(bamChunkReads, observation->readNo);
            RleString *rleString = bamChunkRead->rleRead;
            Symbol symbol = poa->alphabet->convertCharToSymbol(
                    rleString->rleString[observation->offset + observationOffset]);
            int64_t runLength = rleString->repeatCounts[observation->offset + observationOffset];
            bool forward = bamChunkRead->forwardStrand;

            // get correct run length
            runLength -= currentRunLengthIndex * maxRunLength;
            if (runLength < 0) {
                runLength = 0;
            } else if (runLength > maxRunLength) {
                runLength = maxRunLength;
                beforeMaxObservedRunLength = TRUE;
            }

            int64_t pos = PoaFeature_SplitRleWeight_charIndex(maxRunLength, symbol, runLength, forward);
            currFeature->weights[pos] += observation->weight;

        }

        // update currFeature if we're going ot run again
        if (beforeMaxObservedRunLength) {
            currentRunLengthIndex++;
            if (currFeature->nextRunLength != NULL) {
                currFeature = currFeature->nextRunLength;
            } else {
                PoaFeatureSplitRleWeight *prevFeature = currFeature;
                currFeature = PoaFeature_SplitRleWeight_construct(baseFeature->refPosition, baseFeature->insertPosition,
                                                                  currentRunLengthIndex, maxRunLength);
                prevFeature->nextRunLength = currFeature;
                currFeature->weights[PoaFeature_SplitRleWeight_gapIndex(maxRunLength, TRUE)] =
                        baseFeature->weights[PoaFeature_SplitRleWeight_gapIndex(maxRunLength, TRUE)];
                currFeature->weights[PoaFeature_SplitRleWeight_gapIndex(maxRunLength, FALSE)] =
                        baseFeature->weights[PoaFeature_SplitRleWeight_gapIndex(maxRunLength, FALSE)];
            }
        }
    }
}


stList *PoaFeature_getSplitRleWeightFeatures(Poa *poa, stList *bamChunkReads, int64_t maxRunLength) {
    // initialize feature list
    stList *featureList = stList_construct3(0, (void (*)(void *)) PoaFeature_SplitRleWeight_destruct);
    for (int64_t i = 1; i < stList_length(poa->nodes); i++) {
        stList_append(featureList, PoaFeature_SplitRleWeight_construct(i - 1, 0, 0, maxRunLength));
    }

    // for logging (of errors)
    char *logIdentifier = getLogIdentifier();

    // iterate over all positions
    for (int64_t i = 0; i < stList_length(featureList); i++) {

        // get feature and node
        PoaFeatureSplitRleWeight *feature = stList_get(featureList, i);
        PoaNode *node = stList_get(poa->nodes, i +
                                               1); //skip the first poa node, as it's always an 'N', so featureIdx and poaIdx are off by one

        // save run length nodes
        poa_addSplitRunLengthFeaturesForObservations(poa, feature, node->observations, bamChunkReads,
                                                     maxRunLength, 0);

        // Deletes
        if (stList_length(node->deletes) > 0) {

            // iterate over all deletes
            for (int64_t d = 0; d < stList_length(node->deletes); d++) {
                PoaDelete *delete = stList_get(node->deletes, d);

                // Deletes start AFTER the current position, need to add counts/weights to nodes after the current node
                for (int64_t k = 1; k < delete->length; k++) {
                    if (i + k >= stList_length(poa->nodes)) {
                        st_logInfo(" %s Encountered DELETE that occurs after end of POA!\n", logIdentifier);
                        break;
                    }
                    PoaFeatureSplitRleWeight *delFeature = stList_get(featureList, i + k);

                    delFeature->weights[PoaFeature_SplitRleWeight_gapIndex(maxRunLength,
                                                                           TRUE)] += delete->weightForwardStrand;
                    delFeature->weights[PoaFeature_SplitRleWeight_gapIndex(maxRunLength,
                                                                           FALSE)] += delete->weightReverseStrand;
                }
            }
        }

        // Inserts
        if (stList_length(node->inserts) > 0) {

            // iterate over all inserts
            for (int64_t n = 0; n < stList_length(node->inserts); n++) {
                PoaInsert *insert = stList_get(node->inserts, n);

                // handle each insert base
                PoaFeatureSplitRleWeight *prevFeature = feature;
                for (int64_t o = 0; o < strlen(insert->insert->rleString); o++) {

                    // get feature iterator
                    PoaFeatureSplitRleWeight *currFeature = prevFeature->nextInsert;
                    if (currFeature == NULL) {
                        currFeature = PoaFeature_SplitRleWeight_construct(i, o + 1, 0, maxRunLength);
                        prevFeature->nextInsert = currFeature;
                    }

                    // save insert run lengths
                    poa_addSplitRunLengthFeaturesForObservations(poa, currFeature, insert->observations, bamChunkReads,
                                                                 maxRunLength, o);
                }
            }
        }
    }

    free(logIdentifier);
    return featureList;
}

void poa_addChannelRunLengthFeaturesForObservations(Poa *poa, PoaFeatureChannelRleWeight *baseFeature, stList *observations,
                                                    stList *bamChunkReads, const int64_t maxRunLength,
                                                    int64_t observationOffset) {

    PoaFeatureChannelRleWeight *currFeature = baseFeature;
    int64_t currentRunLengthIndex = 0;
    bool beforeMaxObservedRunLength = TRUE;

    while (beforeMaxObservedRunLength) {
        beforeMaxObservedRunLength = FALSE;

        // examine each observation
        for (int64_t i = 0; i < stList_length(observations); i++) {

            // get data
            PoaBaseObservation *observation = stList_get(observations, i);
            BamChunkRead *bamChunkRead = stList_get(bamChunkReads, observation->readNo);
            RleString *rleString = bamChunkRead->rleRead;
            Symbol symbol = poa->alphabet->convertCharToSymbol(
                    rleString->rleString[observation->offset + observationOffset]);
            int64_t runLength = rleString->repeatCounts[observation->offset + observationOffset];
            bool forward = bamChunkRead->forwardStrand;

            // get correct run length
            runLength -= currentRunLengthIndex * maxRunLength;
            if (runLength < 0) {
                runLength = 0;
            } else if (runLength > maxRunLength) {
                runLength = maxRunLength;
                beforeMaxObservedRunLength = TRUE;
            }

            // save weight for both nucleotide totals and runLenght totals
            int64_t nuclPos = PoaFeature_ChannelRleWeight_charNuclIndex(symbol, forward);
            int64_t runLengthPos = PoaFeature_ChannelRleWeight_charRLIndex(maxRunLength, symbol, runLength, forward);
            currFeature->nucleotideWeights[nuclPos] += observation->weight;
            currFeature->runLengthWeights[runLengthPos] += observation->weight;

        }

        // update currFeature if we're going ot run again
        if (beforeMaxObservedRunLength) {
            currentRunLengthIndex++;
            if (currFeature->nextRunLength != NULL) {
                currFeature = currFeature->nextRunLength;
            } else {
                PoaFeatureChannelRleWeight *prevFeature = currFeature;
                currFeature = PoaFeature_ChannelRleWeight_construct(baseFeature->refPosition,
                                                                    baseFeature->insertPosition,
                                                                    currentRunLengthIndex, maxRunLength);
                prevFeature->nextRunLength = currFeature;
                currFeature->nucleotideWeights[PoaFeature_ChannelRleWeight_gapNuclIndex(TRUE)] =
                        baseFeature->nucleotideWeights[PoaFeature_ChannelRleWeight_gapNuclIndex(TRUE)];
                currFeature->nucleotideWeights[PoaFeature_ChannelRleWeight_gapNuclIndex(FALSE)] =
                        baseFeature->nucleotideWeights[PoaFeature_ChannelRleWeight_gapNuclIndex(FALSE)];
            }
        }
    }
}

stList *PoaFeature_getChannelRleWeightFeatures(Poa *poa, stList *bamChunkReads, int64_t maxRunLength) {
    // initialize feature list
    stList *featureList = stList_construct3(0, (void (*)(void *)) PoaFeature_ChannelRleWeight_destruct);
    for (int64_t i = 1; i < stList_length(poa->nodes); i++) {
        stList_append(featureList, PoaFeature_ChannelRleWeight_construct(i - 1, 0, 0, maxRunLength));
    }

    // for logging (of errors)
    char *logIdentifier = getLogIdentifier();

    // iterate over all positions
    for (int64_t i = 0; i < stList_length(featureList); i++) {

        // get feature and node
        PoaFeatureChannelRleWeight *feature = stList_get(featureList, i);
        PoaNode *node = stList_get(poa->nodes, i +
                                               1); //skip the first poa node, as it's always an 'N', so featureIdx and poaIdx are off by one

        // save run length nodes
        poa_addChannelRunLengthFeaturesForObservations(poa, feature, node->observations, bamChunkReads,
                                                       maxRunLength, 0);

        // Deletes
        if (stList_length(node->deletes) > 0) {

            // iterate over all deletes
            for (int64_t d = 0; d < stList_length(node->deletes); d++) {
                PoaDelete *delete = stList_get(node->deletes, d);

                // Deletes start AFTER the current position, need to add counts/weights to nodes after the current node
                for (int64_t k = 1; k < delete->length; k++) {
                    if (i + k >= stList_length(poa->nodes)) {
                        st_logCritical(" %s Encountered DELETE that occurs after end of POA!\n", logIdentifier);
                        break;
                    }
                    PoaFeatureChannelRleWeight *delFeature = stList_get(featureList, i + k);

                    delFeature->nucleotideWeights[PoaFeature_ChannelRleWeight_gapNuclIndex(TRUE)] +=
                            delete->weightForwardStrand;
                    delFeature->nucleotideWeights[PoaFeature_ChannelRleWeight_gapNuclIndex(FALSE)] +=
                            delete->weightReverseStrand;
                }
            }
        }

        // Inserts
        if (stList_length(node->inserts) > 0) {

            // iterate over all inserts
            for (int64_t n = 0; n < stList_length(node->inserts); n++) {
                PoaInsert *insert = stList_get(node->inserts, n);

                // handle each insert base
                PoaFeatureChannelRleWeight *prevFeature = feature;
                for (int64_t o = 0; o < strlen(insert->insert->rleString); o++) {

                    // get feature iterator
                    PoaFeatureChannelRleWeight *currFeature = prevFeature->nextInsert;
                    if (currFeature == NULL) {
                        currFeature = PoaFeature_ChannelRleWeight_construct(i, o + 1, 0, maxRunLength);
                        prevFeature->nextInsert = currFeature;
                    }

                    // save insert run lengths
                    poa_addChannelRunLengthFeaturesForObservations(poa, currFeature, insert->observations,
                                                                   bamChunkReads,
                                                                   maxRunLength, o);
                }
            }
        }
    }

    free(logIdentifier);
    return featureList;
}


void poa_addDiploidRunLengthFeaturesForObservations(Poa *poa, PoaFeatureDiploidRleWeight *baseFeature, stList *observations,
                                                  stList *bamChunkReads, stSet *onHapReads, const int64_t maxRunLength,
                                                  int64_t observationOffset) {


    PoaFeatureDiploidRleWeight *currFeature = baseFeature;
    int64_t currentRunLengthIndex = 0;
    bool beforeMaxObservedRunLength = TRUE;

    while (beforeMaxObservedRunLength) {
        beforeMaxObservedRunLength = FALSE;

        // examine each observation
        for (int64_t i = 0; i < stList_length(observations); i++) {

            // get data
            PoaBaseObservation *observation = stList_get(observations, i);
            BamChunkRead *bamChunkRead = stList_get(bamChunkReads, observation->readNo);
            RleString *rleString = bamChunkRead->rleRead;
            Symbol symbol = poa->alphabet->convertCharToSymbol(rleString->rleString[observation->offset + observationOffset]);
            int64_t runLength = rleString->repeatCounts[observation->offset + observationOffset];
            bool forward = bamChunkRead->forwardStrand;

            // get correct run length
            runLength -= currentRunLengthIndex * maxRunLength;
            if (runLength < 0) {
                runLength = 0;
            } else if (runLength > maxRunLength) {
                runLength = maxRunLength;
                beforeMaxObservedRunLength = TRUE;
            }

            int64_t pos = PoaFeature_DiploidRleWeight_charIndex(maxRunLength, symbol, runLength, forward);
            if (stSet_search(onHapReads, bamChunkRead->readName) != NULL) {
                currFeature->weightsHOn[pos] += observation->weight;
            } else {
                currFeature->weightsHOff[pos] += observation->weight;
            }
        }

        // update currFeature if we're going ot run again
        if (beforeMaxObservedRunLength) {
            currentRunLengthIndex++;
            if (currFeature->nextRunLength != NULL) {
                currFeature = currFeature->nextRunLength;
            } else {
                PoaFeatureDiploidRleWeight *prevFeature = currFeature;
                currFeature = PoaFeature_DiploidRleWeight_construct(baseFeature->refPosition, baseFeature->insertPosition,
                                                                  currentRunLengthIndex, maxRunLength);
                prevFeature->nextRunLength = currFeature;

                currFeature->weightsHOn[PoaFeature_DiploidRleWeight_gapIndex(maxRunLength, TRUE)] =
                        baseFeature->weightsHOn[PoaFeature_DiploidRleWeight_gapIndex(maxRunLength, TRUE)];
                currFeature->weightsHOn[PoaFeature_DiploidRleWeight_gapIndex(maxRunLength, FALSE)] =
                        baseFeature->weightsHOn[PoaFeature_DiploidRleWeight_gapIndex(maxRunLength, FALSE)];

                currFeature->weightsHOff[PoaFeature_DiploidRleWeight_gapIndex(maxRunLength, TRUE)] =
                        baseFeature->weightsHOff[PoaFeature_DiploidRleWeight_gapIndex(maxRunLength, TRUE)];
                currFeature->weightsHOff[PoaFeature_DiploidRleWeight_gapIndex(maxRunLength, FALSE)] =
                        baseFeature->weightsHOff[PoaFeature_DiploidRleWeight_gapIndex(maxRunLength, FALSE)];
            }
        }
    }
}



stList *PoaFeature_getDiploidRleWeightFeatures(Poa *poa, stList *bamChunkReads, stSet *onHapReads,
                                               const int64_t maxRunLength) {
    // initialize feature list
    stList *featureList = stList_construct3(0, (void (*)(void *)) PoaFeature_DiploidRleWeight_destruct);
    for(int64_t i=1; i<stList_length(poa->nodes); i++) {
        stList_append(featureList, PoaFeature_DiploidRleWeight_construct(i - 1, 0, 0, maxRunLength));
    }

    // for logging (of errors)
    char *logIdentifier = getLogIdentifier();

    // iterate over all positions
    for(int64_t i=0; i<stList_length(featureList); i++) {

        // get feature and node
        PoaFeatureDiploidRleWeight* feature = stList_get(featureList, i);
        PoaNode *node = stList_get(poa->nodes, i + 1); //skip the first poa node, as it's always an 'N', so featureIdx and poaIdx are off by one

        // save run length nodes
        poa_addDiploidRunLengthFeaturesForObservations(poa, feature, node->observations, bamChunkReads, onHapReads,
                                                       maxRunLength, 0);

        // Deletes
        if (stList_length(node->deletes) > 0) {

            // iterate over all deletes
            for (int64_t d = 0; d < stList_length(node->deletes); d++) {
                PoaDelete *delete = stList_get(node->deletes, d);

                // Deletes start AFTER the current position, need to add counts/weights to nodes after the current node
                for (int64_t k = 1; k < delete->length; k++) {
                    if (i + k >= stList_length(poa->nodes)) {
                        st_logInfo(" %s Encountered DELETE that occurs after end of POA!\n", logIdentifier);
                        break;
                    }
                    PoaFeatureDiploidRleWeight *delFeature = stList_get(featureList, i + k);
                    stList *observations = delete->observations;

                    // examine each observation for hap and direction
                    for (int64_t o = 0; o < stList_length(observations); o++) {

                        // get observation data
                        PoaBaseObservation *observation = stList_get(observations, o);
                        BamChunkRead *bamChunkRead = stList_get(bamChunkReads, observation->readNo);
                        bool forward = bamChunkRead->forwardStrand;

                        int64_t pos = PoaFeature_DiploidRleWeight_gapIndex(maxRunLength, forward);
                        if (stSet_search(onHapReads, bamChunkRead->readName) != NULL) {
                            delFeature->weightsHOn[pos] += observation->weight;
                        } else {
                            delFeature->weightsHOff[pos] += observation->weight;
                        }
                    }
                }
            }
        }

        // Inserts
        if (stList_length(node->inserts) > 0) {

            // iterate over all inserts
            for (int64_t n = 0; n < stList_length(node->inserts); n++) {
                PoaInsert *insert = stList_get(node->inserts, n);

                // handle each insert base
                PoaFeatureDiploidRleWeight *prevFeature = feature;
                for (int64_t o = 0; o < strlen(insert->insert->rleString); o++) {

                    // get feature iterator
                    PoaFeatureDiploidRleWeight *currFeature = prevFeature->nextInsert;
                    if (currFeature == NULL) {
                        currFeature = PoaFeature_DiploidRleWeight_construct(i, o + 1, 0, maxRunLength);
                        prevFeature->nextInsert = currFeature;
                    }

                    // save insert run lengths
                    poa_addDiploidRunLengthFeaturesForObservations(poa, currFeature, insert->observations,
                            bamChunkReads, onHapReads, maxRunLength, o);
                }
            }
        }
    }

    free(logIdentifier);
    return featureList;
}

void printMEAAlignment2(RleString *X, RleString *Y, stList *alignedPairs) {
    printMEAAlignment(X->rleString, Y->rleString, X->length, Y->length, alignedPairs, X->repeatCounts, Y->repeatCounts);
    fprintf(stderr, "          AlignmentIdentity: %f\n", calculateAlignIdentity(X, Y, alignedPairs));
}

void printMEAAlignment(char *X, char *Y, int64_t lX, int64_t lY, stList *alignedPairs, uint64_t *Xrl, uint64_t *Yrl) {
    if (stList_length(alignedPairs) == 0) {
        return;
    }

    // should we do run lengths
    bool handleRunLength = Xrl != NULL && Yrl != NULL;

    // build strings to print
    int64_t bufferLen = lX + lY;
    char *alnXStr = st_calloc(bufferLen, sizeof(char));
    char *alnYStr = st_calloc(bufferLen, sizeof(char));
    char *alnDesc = st_calloc(bufferLen, sizeof(char));
    char *rlXStr = NULL;
    char *rlYStr = NULL;
    if (handleRunLength) {
        rlXStr = st_calloc(bufferLen, sizeof(char));
        rlYStr = st_calloc(bufferLen, sizeof(char));
    }

    // stats to track
    int64_t nuclMatches = 0;
    int64_t nuclMismatches = 0;
    int64_t nuclXInserts = 0;
    int64_t nuclYInserts = 0;
    int64_t rlMatches = 0;
    int64_t rlMismatches = 0;


    // iterate over alignment
    stListIterator *alignmentItor = stList_getIterator(alignedPairs);
    stIntTuple *prevAlign = NULL;
    stIntTuple *currAlign = stList_getNext(alignmentItor);
    int64_t posX = stIntTuple_get(currAlign, 0);
    int64_t posY = stIntTuple_get(currAlign, 1);
    int64_t origPosX = posX;
    int64_t origPosY = posY;
    int64_t outStrPos = 0;

    while (TRUE) {
        if (currAlign == NULL) break;
        int64_t currAlignPosX = stIntTuple_get(currAlign, 0);
        int64_t currAlignPosY = stIntTuple_get(currAlign, 1);
        // Y gap / X insert
        if (posX < currAlignPosX) {
            alnXStr[outStrPos] = X[posX];
            alnYStr[outStrPos] = '_';
            alnDesc[outStrPos] = ' ';
            if (handleRunLength) {
                rlXStr[outStrPos] = (char) ('0' + Xrl[posX]);
                rlYStr[outStrPos] = ' ';
            }
            posX++;
            nuclXInserts++;
        }

            // X gap / Y insert
        else if (posY < currAlignPosY) {
            alnXStr[outStrPos] = '_';
            alnYStr[outStrPos] = Y[posY];
            alnDesc[outStrPos] = ' ';
            if (handleRunLength) {
                rlXStr[outStrPos] = ' ';
                rlYStr[outStrPos] = (char) ('0' + Yrl[posY]);
            }
            posY++;
            nuclYInserts++;
        }

            // match
        else if (posX == currAlignPosX && posY == currAlignPosY) {
            alnXStr[outStrPos] = X[posX];
            alnYStr[outStrPos] = Y[posY];
            if (handleRunLength) {
                rlXStr[outStrPos] = (char) ('0' + Xrl[posX]);
                rlYStr[outStrPos] = (char) ('0' + Yrl[posY]);
            }
            if (tolower(X[posX]) == tolower(Y[posY])) {
                nuclMatches++;
                alnDesc[outStrPos] = '|';
                if (handleRunLength) {
                    if (Xrl[posX] == Yrl[posY]) {
                        rlMatches++;
                    } else {
                        rlMismatches++;
                        alnDesc[outStrPos] = ':';
                    }
                }
            } else {
                alnDesc[outStrPos] = ' ';
                nuclMismatches++;
            }
            posX++;
            posY++;
            prevAlign = currAlign;
            currAlign = stList_getNext(alignmentItor);
        }

            // should never happen
        else {
            assert(FALSE);
        }

        outStrPos++;
    }

    // print
//    fprintf(stderr, "\n");
    if (handleRunLength) fprintf(stderr, "          %s\n", rlXStr);
    fprintf(stderr, "  %7"PRId64" %s %7"PRId64"\n", origPosX, alnXStr, posX);
    fprintf(stderr, "          %s\n", alnDesc);
    fprintf(stderr, "  %7"PRId64" %s %7"PRId64"\n", origPosY, alnYStr, posY);
    if (handleRunLength) fprintf(stderr, "          %s\n", rlYStr);
    fprintf(stderr, "          Matches:    %"PRId64"\n", nuclMatches);
    if (handleRunLength) {
        fprintf(stderr, "            RL Match: %"PRId64"\n", rlMatches);
        fprintf(stderr, "            RL Miss:  %"PRId64"\n", rlMismatches);
    }
    fprintf(stderr, "          Mismatches: %"PRId64"\n", nuclMismatches);
    fprintf(stderr, "          X Inserts:  %"PRId64"\n", nuclXInserts);
    fprintf(stderr, "          Y Inserts:  %"PRId64"\n", nuclYInserts);
//    fprintf(stderr, "\n");

    // cleanup
    free(alnXStr);
    free(alnYStr);
    free(alnDesc);
    if (handleRunLength) {
        free(rlXStr);
        free(rlYStr);
    }
    stList_destructIterator(alignmentItor);
}


void annotateHelenFeaturesWithTruth(stList *features, HelenFeatureType featureType, stList *trueRefAlignment,
                                    RleString *trueRefRleString, int64_t *firstMatchedFeaure,
                                    int64_t *lastMatchedFeature) {
    /*
     Each index in features represents a character in the final consensus string (which in turn was a single node in the
     final POA which was used to generate the consensus).  This consensus was aligned against the true reference (which
     is reflected in the trueRefRleString).  Items in the true refAlignment are stIntTuples with index 0 being the
     alignment weight (discarded), index 1 being the position in the consensus string (and features), and index 2 being
     the position in the trueRefRleString.  So we can iterate over features and the true alignment to assign truth
     labels to each feature.
     */
    static int FEATURE_POS = 0;
    static int REFERENCE_POS = 1;
    *firstMatchedFeaure = -1;
    *lastMatchedFeature = -1;
    char *logIdentifier = getLogIdentifier();

    // iterate over true ref alignment
    stListIterator *trueRefAlignItor = stList_getIterator(trueRefAlignment);
    stIntTuple *currTrueRefAlign = stList_getNext(trueRefAlignItor);

    // iterate over features
    int64_t trueRefPos = stIntTuple_get(currTrueRefAlign, REFERENCE_POS);
    for (int64_t featureRefPos = 0; featureRefPos < stList_length(features); featureRefPos++) {
        void *feature = stList_get(features, featureRefPos);
        int64_t trueRunLength = -1;
        PoaFeatureSplitRleWeight* srlFeature = NULL;
        PoaFeatureChannelRleWeight* crlFeature = NULL;
        PoaFeatureDiploidRleWeight* drlFeature = NULL;

        int64_t featureInsPos = 0;
        while (feature != NULL) {

            // no more ref bases, everything is gaps
            if (currTrueRefAlign == NULL) {
                switch (featureType) {
                    case HFEAT_SIMPLE_WEIGHT:
                        ((PoaFeatureSimpleWeight *) feature)->label = '_';
                        feature = ((PoaFeatureSimpleWeight *) feature)->nextInsert;
                        break;
                    case HFEAT_SPLIT_RLE_WEIGHT:
                        srlFeature = ((PoaFeatureSplitRleWeight *) feature);
                        while (srlFeature != NULL) {
                            srlFeature->labelChar = '_';
                            srlFeature->labelRunLength = 0;
                            srlFeature = srlFeature->nextRunLength;
                        }
                        feature = ((PoaFeatureSplitRleWeight *) feature)->nextInsert;
                        break;
                    case HFEAT_CHANNEL_RLE_WEIGHT:
                        crlFeature = ((PoaFeatureChannelRleWeight *) feature);
                        while (crlFeature != NULL) {
                            crlFeature->labelChar = '_';
                            crlFeature->labelRunLength = 0;
                            crlFeature = crlFeature->nextRunLength;
                        }
                        feature = ((PoaFeatureChannelRleWeight *) feature)->nextInsert;
                        break;
                    case HFEAT_DIPLOID_RLE_WEIGHT:
                        drlFeature = ((PoaFeatureDiploidRleWeight*)feature);
                        while (drlFeature != NULL) {
                            drlFeature->labelChar = '_';
                            drlFeature->labelRunLength = 0;
                            drlFeature = drlFeature->nextRunLength;
                        }
                        feature = ((PoaFeatureDiploidRleWeight*)feature)->nextInsert;
                        break;
                    default:
                        st_errAbort("Unhandled FeatureType!\n");
                }
                continue;
            }

            // sanity checks
            assert(stIntTuple_get(currTrueRefAlign, FEATURE_POS) >= featureRefPos &&
                   stIntTuple_get(currTrueRefAlign, REFERENCE_POS) >= trueRefPos);

            // match
            if (stIntTuple_get(currTrueRefAlign, FEATURE_POS) == featureRefPos &&
                stIntTuple_get(currTrueRefAlign, REFERENCE_POS) == trueRefPos) {
                st_logDebug(
                        " %s LABEL MATCH  %c trueRefPos:%"PRId64" featureRefPos:%"PRId64" featureInsPos:%"PRId64"\n",
                        logIdentifier, featureInsPos == 0 ? ' ' : 'I', trueRefPos, featureRefPos, featureInsPos);
                // save label (based on feature type)
                switch (featureType) {
                    case HFEAT_SIMPLE_WEIGHT:
                        ((PoaFeatureSimpleWeight *) feature)->label = trueRefRleString->rleString[trueRefPos];
                        break;
                    case HFEAT_SPLIT_RLE_WEIGHT:
                        srlFeature = ((PoaFeatureSplitRleWeight *) feature);
                        trueRunLength = trueRefRleString->repeatCounts[trueRefPos];
                        while (srlFeature != NULL) {
                            srlFeature->labelChar = trueRefRleString->rleString[trueRefPos];
                            if (trueRunLength <= 0) {
                                srlFeature->labelRunLength = 0;
                            } else if (trueRunLength > srlFeature->maxRunLength) {
                                srlFeature->labelRunLength = srlFeature->maxRunLength;
                            } else {
                                srlFeature->labelRunLength = trueRunLength;
                            }
                            trueRunLength -= srlFeature->maxRunLength;
                            srlFeature = srlFeature->nextRunLength;
                        }
                        break;
                    case HFEAT_CHANNEL_RLE_WEIGHT:
                        crlFeature = ((PoaFeatureChannelRleWeight *) feature);
                        trueRunLength = trueRefRleString->repeatCounts[trueRefPos];
                        while (crlFeature != NULL) {
                            crlFeature->labelChar = trueRefRleString->rleString[trueRefPos];
                            if (trueRunLength <= 0) {
                                crlFeature->labelRunLength = 0;
                            } else if (trueRunLength > crlFeature->maxRunLength) {
                                crlFeature->labelRunLength = crlFeature->maxRunLength;
                            } else {
                                crlFeature->labelRunLength = trueRunLength;
                            }
                            trueRunLength -= crlFeature->maxRunLength;
                            crlFeature = crlFeature->nextRunLength;
                        }
                        break;
                    case HFEAT_DIPLOID_RLE_WEIGHT:
                        drlFeature = ((PoaFeatureDiploidRleWeight*)feature);
                        trueRunLength = trueRefRleString->repeatCounts[trueRefPos];
                        while (drlFeature != NULL) {
                            drlFeature->labelChar = trueRefRleString->rleString[trueRefPos];
                            if (trueRunLength <= 0) {
                                drlFeature->labelRunLength = 0;
                            } else if (trueRunLength > drlFeature->maxRunLength) {
                                drlFeature->labelRunLength = drlFeature->maxRunLength;
                            } else {
                                drlFeature->labelRunLength = trueRunLength;
                            }
                            trueRunLength -= drlFeature->maxRunLength;
                            drlFeature = drlFeature->nextRunLength;
                        }
                        break;
                    default:
                        st_errAbort("Unhandled FeatureType!\n");
                }

                // iterate
                trueRefPos++;
                currTrueRefAlign = stList_getNext(trueRefAlignItor);
                // handle first and last match
                if (featureInsPos == 0) {
                    if (*firstMatchedFeaure == -1) {
                        *firstMatchedFeaure = featureRefPos;
                    }
                    *lastMatchedFeature = featureRefPos;
                }
            }

                // insert
            else if (trueRefPos < stIntTuple_get(currTrueRefAlign, REFERENCE_POS)) {
                st_logDebug(
                        " %s LABEL INSERT %c trueRefPos:%"PRId64" featureRefPos:%"PRId64" featureInsPos:%"PRId64"\n",
                        logIdentifier, featureInsPos == 0 ? ' ' : 'I', trueRefPos, featureRefPos, featureInsPos);
                // apply label
                switch (featureType) {
                    case HFEAT_SIMPLE_WEIGHT:
                        ((PoaFeatureSimpleWeight *) feature)->label = trueRefRleString->rleString[trueRefPos];
                        break;
                    case HFEAT_SPLIT_RLE_WEIGHT:
                        srlFeature = ((PoaFeatureSplitRleWeight *) feature);
                        trueRunLength = trueRefRleString->repeatCounts[trueRefPos];
                        while (srlFeature != NULL) {
                            srlFeature->labelChar = trueRefRleString->rleString[trueRefPos];
                            if (trueRunLength <= 0) {
                                srlFeature->labelRunLength = 0;
                            } else if (trueRunLength > srlFeature->maxRunLength) {
                                srlFeature->labelRunLength = srlFeature->maxRunLength;
                            } else {
                                srlFeature->labelRunLength = trueRunLength;
                            }
                            trueRunLength -= srlFeature->maxRunLength;
                            srlFeature = srlFeature->nextRunLength;
                        }
                        break;
                    case HFEAT_CHANNEL_RLE_WEIGHT:
                        crlFeature = ((PoaFeatureChannelRleWeight *) feature);
                        trueRunLength = trueRefRleString->repeatCounts[trueRefPos];
                        while (crlFeature != NULL) {
                            crlFeature->labelChar = trueRefRleString->rleString[trueRefPos];
                            if (trueRunLength <= 0) {
                                crlFeature->labelRunLength = 0;
                            } else if (trueRunLength > crlFeature->maxRunLength) {
                                crlFeature->labelRunLength = crlFeature->maxRunLength;
                            } else {
                                crlFeature->labelRunLength = trueRunLength;
                            }
                            trueRunLength -= crlFeature->maxRunLength;
                            crlFeature = crlFeature->nextRunLength;
                        }
                        break;
                    case HFEAT_DIPLOID_RLE_WEIGHT:
                        drlFeature = ((PoaFeatureDiploidRleWeight*)feature);
                        trueRunLength = trueRefRleString->repeatCounts[trueRefPos];
                        while (drlFeature != NULL) {
                            drlFeature->labelChar = trueRefRleString->rleString[trueRefPos];
                            if (trueRunLength <= 0) {
                                drlFeature->labelRunLength = 0;
                            } else if (trueRunLength > drlFeature->maxRunLength) {
                                drlFeature->labelRunLength = drlFeature->maxRunLength;
                            } else {
                                drlFeature->labelRunLength = trueRunLength;
                            }
                            trueRunLength -= drlFeature->maxRunLength;
                            drlFeature = drlFeature->nextRunLength;
                        }
                        break;
                    default:
                        st_errAbort("Unhandled FeatureType!\n");
                }
                trueRefPos++;
            }

                // delete
            else if (featureRefPos < stIntTuple_get(currTrueRefAlign, FEATURE_POS)) {
                st_logDebug(
                        " %s LABEL DELETE %c trueRefPos:%"PRId64" featureRefPos:%"PRId64" featureInsPos:%"PRId64"\n",
                        logIdentifier, featureInsPos == 0 ? ' ' : 'I', trueRefPos, featureRefPos, featureInsPos);
                // apply label
                switch (featureType) {
                    case HFEAT_SIMPLE_WEIGHT:
                        ((PoaFeatureSimpleWeight *) feature)->label = '_';
                        break;
                    case HFEAT_SPLIT_RLE_WEIGHT:
                        srlFeature = ((PoaFeatureSplitRleWeight *) feature);
                        while (srlFeature != NULL) {
                            srlFeature->labelChar = '_';
                            srlFeature->labelRunLength = 0;
                            srlFeature = srlFeature->nextRunLength;
                        }
                        break;
                    case HFEAT_CHANNEL_RLE_WEIGHT:
                        crlFeature = ((PoaFeatureChannelRleWeight *) feature);
                        while (crlFeature != NULL) {
                            crlFeature->labelChar = '_';
                            crlFeature->labelRunLength = 0;
                            crlFeature = crlFeature->nextRunLength;
                        }
                        break;
                    case HFEAT_DIPLOID_RLE_WEIGHT:
                        drlFeature = ((PoaFeatureDiploidRleWeight*)feature);
                        while (drlFeature != NULL) {
                            drlFeature->labelChar = '_';
                            drlFeature->labelRunLength = 0;
                            drlFeature = drlFeature->nextRunLength;
                        }
                        break;
                    default:
                        st_errAbort("Unhandled FeatureType!\n");
                }
            }

                // programmer error
            else {
                st_errAbort("Unhandled case annotating features with true reference characters!\n");
            }

            // always iterate over insert features
            switch (featureType) {
                case HFEAT_SIMPLE_WEIGHT:
                    feature = ((PoaFeatureSimpleWeight *) feature)->nextInsert;
                    break;
                case HFEAT_SPLIT_RLE_WEIGHT:
                    feature = ((PoaFeatureSplitRleWeight *) feature)->nextInsert;
                    break;
                case HFEAT_CHANNEL_RLE_WEIGHT:
                    feature = ((PoaFeatureChannelRleWeight *) feature)->nextInsert;
                    break;
                case HFEAT_DIPLOID_RLE_WEIGHT:
                    feature = ((PoaFeatureDiploidRleWeight*)feature)->nextInsert;
                    break;
                default:
                    st_errAbort("Unhandled FeatureType!\n");
            }
            featureInsPos++;
        }

        // this catches any true inserts which are not present in the poa / feature list
        while (currTrueRefAlign != NULL &&
               featureRefPos < stIntTuple_get(currTrueRefAlign, FEATURE_POS) &&
               trueRefPos < stIntTuple_get(currTrueRefAlign, REFERENCE_POS)) {
            trueRefPos++;
        }
    }

    stList_destructIterator(trueRefAlignItor);
    free(logIdentifier);
}

void PoaFeature_writeHelenFeatures(HelenFeatureType type, Poa *poa, stList *bamChunkReads,
                                   char *outputFileBase, BamChunk *bamChunk, stList *trueRefAlignment,
                                   RleString *consensusRleString,
                                   RleString *trueRefRleString, bool fullFeatureOutput, int64_t splitWeightMaxRunLength,
                                   HelenFeatureHDF5FileInfo **helenHDF5Files) {
    // prep
    int64_t firstMatchedFeature = -1;
    int64_t lastMatchedFeature = -1;
    stList *features = NULL;
    bool outputLabels = trueRefAlignment != NULL && trueRefRleString != NULL;

# ifdef _OPENMP
    int64_t threadIdx = omp_get_thread_num();
# else
    int64_t threadIdx = 0;
# endif

    // handle differently based on type
    switch (type) {
        case HFEAT_SIMPLE_WEIGHT :
            // get features
            features = PoaFeature_getSimpleWeightFeatures(poa, bamChunkReads);
            firstMatchedFeature = 0;
            lastMatchedFeature = stList_length(features) - 1;

            // get truth (if we have it)
            if (outputLabels) {
                annotateHelenFeaturesWithTruth(features, type, trueRefAlignment, trueRefRleString,
                                               &firstMatchedFeature, &lastMatchedFeature);
            }

            writeSimpleWeightHelenFeaturesHDF5(poa->alphabet, helenHDF5Files[threadIdx], outputFileBase, bamChunk,
                                               outputLabels, features, firstMatchedFeature, lastMatchedFeature);

            break;

        case HFEAT_SPLIT_RLE_WEIGHT:
            // get features
            features = PoaFeature_getSplitRleWeightFeatures(poa, bamChunkReads, splitWeightMaxRunLength);
            firstMatchedFeature = 0;
            lastMatchedFeature = stList_length(features) - 1;

            // get truth (if we have it)
            if (outputLabels) {
                annotateHelenFeaturesWithTruth(features, type, trueRefAlignment, trueRefRleString,
                                               &firstMatchedFeature, &lastMatchedFeature);
            }

            writeSplitRleWeightHelenFeaturesHDF5(poa->alphabet, helenHDF5Files[threadIdx],
                                                 outputFileBase, bamChunk, outputLabels, features, firstMatchedFeature, lastMatchedFeature,
                                                 splitWeightMaxRunLength);
            break;

        case HFEAT_CHANNEL_RLE_WEIGHT:
            // get features
            features = PoaFeature_getChannelRleWeightFeatures(poa, bamChunkReads, splitWeightMaxRunLength);
            firstMatchedFeature = 0;
            lastMatchedFeature = stList_length(features) - 1;

            // get truth (if we have it)
            if (outputLabels) {
                annotateHelenFeaturesWithTruth(features, type, trueRefAlignment, trueRefRleString,
                                               &firstMatchedFeature, &lastMatchedFeature);
            }

            writeChannelRleWeightHelenFeaturesHDF5(poa->alphabet, helenHDF5Files[threadIdx],
                                                   outputFileBase, bamChunk, outputLabels, features, firstMatchedFeature,
                                                   lastMatchedFeature, splitWeightMaxRunLength);
            break;

        default:
            st_errAbort("Unhandled HELEN feature type!\n");
    }

    //cleanup
    stList_destruct(features);
}

void PoaFeature_writeDiploidHelenFeatures(HelenFeatureType type, stList *bamChunkReads, char *outputFileBase,
                                          BamChunk *bamChunk, Poa *poaH1, Poa *poaH2, stSet *readsInH1,
                                          stSet *readsInH2, stList *trueRefAlignmentsToH1, stList *trueRefAlignmentsToH2,
                                          int64_t maxRunLength, HelenFeatureHDF5FileInfo **helenHDF5Files) {

    // prep
    int64_t firstMatchedFeatureH1 = -1;
    int64_t firstMatchedFeatureH2 = -1;
    int64_t lastMatchedFeatureH1 = -1;
    int64_t lastMatchedFeatureH2 = -1;
    stList *featuresH1 = NULL;
    stList *featuresH2 = NULL;
    bool outputLabels = trueRefAlignmentsToH1 != NULL && trueRefAlignmentsToH2 != NULL &&
            (stList_length(trueRefAlignmentsToH1) > 0 || stList_length(trueRefAlignmentsToH2) > 0);

# ifdef _OPENMP
    int64_t threadIdx = omp_get_thread_num();
# else
    int64_t threadIdx = 0;
# endif

    // handle differently based on type
    switch (type) {
        case HFEAT_DIPLOID_RLE_WEIGHT :
            // get features
            featuresH1 = PoaFeature_getDiploidRleWeightFeatures(poaH1, bamChunkReads, readsInH1, maxRunLength);
            featuresH2 = PoaFeature_getDiploidRleWeightFeatures(poaH2, bamChunkReads, readsInH2, maxRunLength);

            firstMatchedFeatureH1 = 0;
            firstMatchedFeatureH2 = 0;
            lastMatchedFeatureH1 = stList_length(featuresH1) - 1;
            lastMatchedFeatureH2 = stList_length(featuresH2) - 1;

            // get truth (if we have it)
            if (outputLabels) {
                for (int64_t i = 0; i < stList_length(trueRefAlignmentsToH1); i++) {
                    HelenFeatureTruthAlignment *hfta = stList_get(trueRefAlignmentsToH1, i);
                    annotateHelenFeaturesWithTruth(featuresH1, type, hfta->alignedPairs, hfta->truthSequence,
                                                   &firstMatchedFeatureH1, &lastMatchedFeatureH1);

                    writeDiploidRleWeightHelenFeaturesHDF5(poaH1->alphabet, helenHDF5Files[threadIdx], outputFileBase,
                            bamChunk, 1, i, featuresH1, firstMatchedFeatureH1, lastMatchedFeatureH1, maxRunLength);
                }

                for (int64_t i = 0; i < stList_length(trueRefAlignmentsToH2); i++) {
                    HelenFeatureTruthAlignment *hfta = stList_get(trueRefAlignmentsToH2, i);
                    annotateHelenFeaturesWithTruth(featuresH2, type, hfta->alignedPairs, hfta->truthSequence,
                                                   &firstMatchedFeatureH2, &lastMatchedFeatureH2);

                    writeDiploidRleWeightHelenFeaturesHDF5(poaH1->alphabet, helenHDF5Files[threadIdx], outputFileBase,
                            bamChunk, 2, i, featuresH2, firstMatchedFeatureH2, lastMatchedFeatureH2, maxRunLength);
                }
            } else {
                writeDiploidRleWeightHelenFeaturesHDF5(poaH1->alphabet, helenHDF5Files[threadIdx], outputFileBase,
                        bamChunk, 1, -1, featuresH1, firstMatchedFeatureH1, lastMatchedFeatureH1, maxRunLength);
                writeDiploidRleWeightHelenFeaturesHDF5(poaH1->alphabet, helenHDF5Files[threadIdx], outputFileBase,
                        bamChunk, 2, -1, featuresH2, firstMatchedFeatureH2, lastMatchedFeatureH2, maxRunLength);

            }


            break;

        default:
            st_errAbort("Unhandled HELEN feature type for diploid!\n");
    }

    //cleanup
    stList_destruct(featuresH1);
    stList_destruct(featuresH2);
}

HelenFeatureTruthAlignment *HelenFeatureTruthAlignment_construct(int64_t startPosInclusive, int64_t endPosExclusive,
                                                                 stList *alignedPairs, RleString *truthSequence) {
    HelenFeatureTruthAlignment *hfta = st_calloc(1, sizeof(HelenFeatureTruthAlignment));
    hfta->startPosIncl = startPosInclusive;
    hfta->endPosExcl = endPosExclusive;
    hfta->alignedPairs = alignedPairs;
    hfta->truthSequence = truthSequence;
    return hfta;
}
void HelenFeatureTruthAlignment_destruct(HelenFeatureTruthAlignment *hfta) {
    stList_destruct(hfta->alignedPairs);
    rleString_destruct(hfta->truthSequence);
    free(hfta);
}
int HelenFeatureTruthAlignment_cmp(const HelenFeatureTruthAlignment *hfta1, const HelenFeatureTruthAlignment *hfta2) {
    if (hfta1->startPosIncl < hfta2->startPosIncl) {
        return -1;
    } else if (hfta1->startPosIncl > hfta2->startPosIncl) {
        return 1;
    } else {
        return hfta1->endPosExcl < hfta2->endPosExcl ? -1 : 1;
    }
}


stList *alignConsensusAndTruthRLE(RleString *consensusStr, RleString *truthStr, double *score, PolishParams *polishParams) {


    // Symbol strings
    SymbolString sX = rleString_constructSymbolString(consensusStr, 0, consensusStr->length, polishParams->alphabet,
                                                      TRUE, (uint64_t) polishParams->repeatSubMatrix->maximumRepeatLength);
    SymbolString sY = rleString_constructSymbolString(truthStr, 0, truthStr->length, polishParams->alphabet,
                                                      TRUE, (uint64_t) polishParams->repeatSubMatrix->maximumRepeatLength);
    uint16_t apScore = 0;

    // Run the alignment
    stList *alignedPairs = stList_construct3(0, (void(*)(void*))stIntTuple_destruct);
    stList *gapXPairs = stList_construct3(0, (void(*)(void*))stIntTuple_destruct);
    stList *gapYPairs = stList_construct3(0, (void(*)(void*))stIntTuple_destruct);

    getAlignedPairsWithIndels(polishParams->stateMachineForForwardStrandRead, sX, sY,
            polishParams->p, &alignedPairs, &gapXPairs, &gapYPairs, TRUE, TRUE);
    stList *meaAlignedPairs = getMaximalExpectedAccuracyPairwiseAlignment(alignedPairs, gapXPairs, gapYPairs,
                                                                          sX.length, sY.length, score, polishParams->p);

    // refactor
    stList *finalAlignedPairs = stList_construct3(0, (void(*)(void*))stIntTuple_destruct);
    for (int64_t i = 0; i < stList_length(meaAlignedPairs); i++) {
        stIntTuple *ap = stList_get(meaAlignedPairs, i);
        stList_append(finalAlignedPairs, stIntTuple_construct3(stIntTuple_get(ap, 1),
                                                               stIntTuple_get(ap, 2), stIntTuple_get(ap, 0)));
    }

    // Cleanup
    stList_destruct(meaAlignedPairs);
    stList_destruct(alignedPairs);
    stList_destruct(gapXPairs);
    stList_destruct(gapYPairs);
    symbolString_destruct(sX);
    symbolString_destruct(sY);

    return finalAlignedPairs;
}


stList *alignConsensusAndTruthRLEWithKmerAnchors(RleString *consensusStr, RleString *truthStr, double *score,
                                                 PolishParams *polishParams) {
    // Symbol strings
    uint64_t maxRL = polishParams->useRunLengthEncoding ? (uint64_t) polishParams->repeatSubMatrix->maximumRepeatLength : 2;
    SymbolString sX = rleString_constructSymbolString(consensusStr, 0, consensusStr->length, polishParams->alphabet,
                                                      TRUE, maxRL);
    SymbolString sY = rleString_constructSymbolString(truthStr, 0, truthStr->length, polishParams->alphabet,
                                                      TRUE, maxRL);
    uint16_t apScore = 0;

    // Run the alignment
    stList *alignedPairs = NULL;
    stList *gapXPairs = NULL;
    stList *gapYPairs = NULL;
    stList *anchorPairs = getKmerAlignmentAnchors(sX, sY, (uint64_t) polishParams->p->diagonalExpansion);

    // quick fail
    int64_t minLength = (consensusStr->length < truthStr->length ? consensusStr : truthStr)->length;
    double apRatio = 1.0 * stList_length(anchorPairs) / minLength;
    if (apRatio < .2) {
        char *logIdentifer = getLogIdentifier();
        st_logInfo(" %s got %"PRId64" anchor pairs for min seq len %"PRId64" (%f), not attempting alignment.\n",
                logIdentifer, stList_length(anchorPairs), minLength, apRatio);
        free(logIdentifer);
        stList_destruct(anchorPairs);
        return stList_construct();
    }

    time_t apTime = time(NULL);
    getAlignedPairsWithIndelsUsingAnchors(polishParams->stateMachineForForwardStrandRead, sX, sY, anchorPairs,
                                          polishParams->p, &alignedPairs, &gapXPairs, &gapYPairs, FALSE, FALSE);
    time_t meaTime = time(NULL);
    stList *meaAlignedPairs = getMaximalExpectedAccuracyPairwiseAlignment(alignedPairs, gapXPairs, gapYPairs,
                                                                          sX.length, sY.length, score, polishParams->p);
    char *logIdentifer = getLogIdentifier();
    st_logInfo(" %s Sequence alignment (seq len %"PRId64", %.3f K-AP ratio) got %"PRId64" MEA aligned pairs in %ds and MEA in %ds\n",
            logIdentifer, minLength, apRatio, stList_length(meaAlignedPairs), meaTime-apTime, time(NULL)-meaTime);

    // refactor
    stList *finalAlignedPairs = stList_construct3(0, (void(*)(void*))stIntTuple_destruct);
    for (int64_t i = 0; i < stList_length(meaAlignedPairs); i++) {
        stIntTuple *ap = stList_get(meaAlignedPairs, i);
        stList_append(finalAlignedPairs, stIntTuple_construct3(stIntTuple_get(ap, 1),
                                                               stIntTuple_get(ap, 2), stIntTuple_get(ap, 0)));
    }

    // Cleanup
    stList_destruct(meaAlignedPairs);
    stList_destruct(alignedPairs);
    stList_destruct(gapXPairs);
    stList_destruct(gapYPairs);
    stList_destruct(anchorPairs);
    symbolString_destruct(sX);
    symbolString_destruct(sY);
    free(logIdentifer);

    return finalAlignedPairs;
}


stList *alignConsensusAndTruthRLEWithSSWAnchors(RleString *consensusStr, RleString *truthStr, double *score,
                                                 PolishParams *polishParams) {


    // Symbol strings
    SymbolString sX = rleString_constructSymbolString(consensusStr, 0, consensusStr->length, polishParams->alphabet,
                                                      TRUE, (uint64_t) polishParams->repeatSubMatrix->maximumRepeatLength);
    SymbolString sY = rleString_constructSymbolString(truthStr, 0, truthStr->length, polishParams->alphabet,
                                                      TRUE, (uint64_t) polishParams->repeatSubMatrix->maximumRepeatLength);
    uint16_t apScore = 0;

    // Run the alignment
    stList *alignedPairs = stList_construct3(0, (void(*)(void*))stIntTuple_destruct);
    stList *gapXPairs = stList_construct3(0, (void(*)(void*))stIntTuple_destruct);
    stList *gapYPairs = stList_construct3(0, (void(*)(void*))stIntTuple_destruct);
    stList *anchorPairsTmp = alignConsensusAndTruthSSW(consensusStr->rleString, truthStr->rleString, &apScore);
    stList *anchorPairs = stList_construct3(0, (void(*)(void*))stIntTuple_destruct);
    if (stList_length(anchorPairs) >= 2) {
        stIntTuple *first = stList_get(anchorPairsTmp, 0);
        stIntTuple *last = stList_get(anchorPairsTmp, stList_length(anchorPairs) - 1);
        stList_append(anchorPairs, stIntTuple_construct3(stIntTuple_get(first, 0), stIntTuple_get(first, 1), PAIR_ALIGNMENT_PROB_1));
        stList_append(anchorPairs, stIntTuple_construct3(stIntTuple_get(last, 0), stIntTuple_get(last, 1), PAIR_ALIGNMENT_PROB_1));
    }
    /*for (int64_t i = 0; i<stList_length(anchorPairs); i++) {
        int64_t *ap = stList_get(anchorPairs, i);
        ap[3] = PAIR_ALIGNMENT_PROB_1;
    }*/


    getAlignedPairsWithIndelsUsingAnchors(polishParams->stateMachineForForwardStrandRead, sX, sY, anchorPairsTmp,
                                          polishParams->p, &alignedPairs, &gapXPairs, &gapYPairs, TRUE, TRUE);
    stList *meaAlignedPairs = getMaximalExpectedAccuracyPairwiseAlignment(alignedPairs, gapXPairs, gapYPairs,
                                                                          sX.length, sY.length, score, polishParams->p);

    // refactor
    stList *finalAlignedPairs = stList_construct3(0, (void(*)(void*))stIntTuple_destruct);
    for (int64_t i = 0; i < stList_length(meaAlignedPairs); i++) {
        stIntTuple *ap = stList_get(meaAlignedPairs, i);
        stList_append(finalAlignedPairs, stIntTuple_construct3(stIntTuple_get(ap, 1),
                                                               stIntTuple_get(ap, 2), stIntTuple_get(ap, 0)));
    }

    // Cleanup
    stList_destruct(meaAlignedPairs);
    stList_destruct(alignedPairs);
    stList_destruct(gapXPairs);
    stList_destruct(gapYPairs);
    stList_destruct(anchorPairs);
    stList_destruct(anchorPairsTmp);
    symbolString_destruct(sX);
    symbolString_destruct(sY);

    return finalAlignedPairs;
}

stList *alignConsensusAndTruthCPECAN(char *consensusStr, char *truthStr, double *score, PolishParams *polishParams) {


    // Symbol strings
    SymbolString sX = symbolString_construct(consensusStr, 0, strlen(consensusStr), polishParams->alphabet);
    SymbolString sY = symbolString_construct(truthStr, 0, strlen(truthStr), polishParams->alphabet);
    uint16_t apScore = 0;

    // Use default state machine for alignment
    StateMachine *sM = stateMachine3_constructNucleotide(threeState);

    // Run the alignment
    stList *alignedPairs = stList_construct3(0, (void(*)(void*))stIntTuple_destruct);
    stList *gapXPairs = stList_construct3(0, (void(*)(void*))stIntTuple_destruct);
    stList *gapYPairs = stList_construct3(0, (void(*)(void*))stIntTuple_destruct);

    getAlignedPairsWithIndels(sM, sX, sY, polishParams->p, &alignedPairs, &gapXPairs,
            &gapYPairs, TRUE, TRUE);
    stList *meaAlignedPairs = getMaximalExpectedAccuracyPairwiseAlignment(alignedPairs, gapXPairs, gapYPairs,
            sX.length, sY.length, score, polishParams->p);

    // refactor
    stList *finalAlignedPairs = stList_construct3(0, (void(*)(void*))stIntTuple_destruct);
    for (int64_t i = 0; i < stList_length(meaAlignedPairs); i++) {
        stIntTuple *ap = stList_get(meaAlignedPairs, i);
        stList_append(finalAlignedPairs, stIntTuple_construct3(stIntTuple_get(ap, 1),
                stIntTuple_get(ap, 2), stIntTuple_get(ap, 0)));
    }

    // Cleanup
    stList_destruct(meaAlignedPairs);
    stList_destruct(alignedPairs);
    stList_destruct(gapXPairs);
    stList_destruct(gapYPairs);
    symbolString_destruct(sX);
    symbolString_destruct(sY);
    stateMachine_destruct(sM);

    return finalAlignedPairs;
}




// this function taken from https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library/blob/master/src/example.c
stList *alignConsensusAndTruthSSW(char *consensusStr, char *truthStr, uint16_t *score) {

    int64_t l, m, k;
    uint8_t match = 2, mismatch = 2, gap_open = 3, gap_extension = 1;    // default parameters for genome sequence alignment

    int64_t consensusLen = strlen(consensusStr);
    int64_t truthLen = strlen(truthStr);

    s_profile *profile;

    int8_t *num = st_calloc(consensusLen, sizeof(int8_t));    // the read sequence represented in numbers
    int8_t *ref_num = st_calloc(truthLen, sizeof(int8_t));  // the read sequence represented in numbers
    s_align *result;

    /* This table is used to transform nucleotide letters into numbers. */
    static const int8_t nt_table[128] = {
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
            4, 4, 4, 4, 3, 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
            4, 4, 4, 4, 3, 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
    };

    // initialize scoring matrix for genome sequences
    //  A  C  G  T	N (or other ambiguous code)
    //  2 -2 -2 -2 	0	A
    // -2  2 -2 -2 	0	C
    // -2 -2  2 -2 	0	G
    // -2 -2 -2  2 	0	T
    //	0  0  0  0  0	N (or other ambiguous code)
    int8_t *mat = (int8_t *) calloc(25, sizeof(int8_t));
    for (l = k = 0; l < 4; ++l) {
        for (m = 0; m < 4; ++m)
            mat[k++] = (int8_t) (l == m ? match : -mismatch);    /* weight_match : -weight_mismatch */
        mat[k++] = 0; // ambiguous base: no penalty
    }
    for (m = 0; m < 5; ++m) mat[k++] = 0;

    for (m = 0; m < consensusLen; ++m) num[m] = nt_table[(int) consensusStr[m]];
    profile = ssw_init(num, (int32_t) consensusLen, mat, 5, 2);
    for (m = 0; m < truthLen; ++m) ref_num[m] = nt_table[(int) truthStr[m]];

    // Only the 8 bit of the flag is setted. ssw_align will always return the best alignment beginning position and cigar.
    result = ssw_align (profile, ref_num, (int32_t) truthLen, gap_open, gap_extension, 1, 0, 0, 15);
    *score = result->score1;

    // Convert from cigar to aligned pairs
    int32_t consensusPos = result->read_begin1;
    int32_t truthPos = result->ref_begin1;
    stList *alignedPairs = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);
    if (result->cigar) {
        for (int32_t cigIdx = 0; cigIdx < result->cigarLen; cigIdx++) {
            char letter = cigar_int_to_op(result->cigar[cigIdx]);
            uint32_t length = cigar_int_to_len(result->cigar[cigIdx]);

            if (letter == 'M') {
                for (int32_t matchIdx = 0; matchIdx < length; matchIdx++) {
                    stList_append(alignedPairs, stIntTuple_construct3(consensusPos, truthPos, 0));
                    consensusPos++;
                    truthPos++;
                }
            } else if (letter == 'I') {
                consensusPos += length;
            } else if (letter == 'D') {
                truthPos += length;
            } else {
                st_errAbort("Got unexpected alignment character '%c' aligning truth and consensus", letter);
            }
        }
    }

    // cleanup
    align_destroy(result);
    init_destroy(profile);
    free(mat);
    free(ref_num);
    free(num);

    return alignedPairs;
}

#define HDF5_FEATURE_SIZE 1000

double **getTwoDArrayDouble(int64_t rowCount, int64_t columnCount, bool zeroValues) {
    double **array = st_calloc(rowCount, sizeof(double *));
    array[0] = (double *) st_calloc(columnCount * rowCount, sizeof(double));
    if (zeroValues) {
        for (int64_t i = 0; i < columnCount * rowCount; i++) {
            (*array)[i] = 0.0;
        }
    }
    for (int64_t i = 1; i < rowCount; i++) {
        array[i] = array[0] + i * columnCount;
    }
    return array;
}

float **getTwoDArrayFloat(int64_t rowCount, int64_t columnCount) {
    float **array = st_calloc(rowCount, sizeof(float *));
    array[0] = (float *) st_calloc(columnCount * rowCount, sizeof(float));
    for (int64_t i = 1; i < rowCount; i++) {
        array[i] = array[0] + i * columnCount;
    }
    return array;
}

uint32_t **getTwoDArrayUInt32(int64_t rowCount, int64_t columnCount) {
    uint32_t **array = st_calloc(rowCount, sizeof(uint32_t *));
    array[0] = (uint32_t *) st_calloc(columnCount * rowCount, sizeof(uint32_t));
    for (int64_t i = 1; i < rowCount; i++) {
        array[i] = array[0] + i * columnCount;
    }
    return array;
}

uint8_t **getTwoDArrayUInt8(int64_t rowCount, int64_t columnCount) {
    uint8_t **array = st_calloc(rowCount, sizeof(uint8_t *));
    array[0] = (uint8_t *) st_calloc(columnCount * rowCount, sizeof(uint8_t));
    for (int64_t i = 1; i < rowCount; i++) {
        array[i] = array[0] + i * columnCount;
    }
    return array;
}

uint8_t ***getThreeDArrayUInt8(int64_t depthCount, int64_t rowCount, int64_t columnCount) {
    uint8_t ***array = st_calloc(depthCount, sizeof(uint8_t **));
    array[0] = (uint8_t **) st_calloc(depthCount * rowCount, sizeof(uint8_t *));
    array[0][0] = (uint8_t *) st_calloc(depthCount * rowCount * columnCount, sizeof(uint8_t));
    for (int64_t i = 0; i < depthCount; i++) {
        array[i] = array[0] + i * rowCount;
        for (int64_t j = 0; j < rowCount; j++) {
            array[i][j] = (array[0][0] + i * rowCount * columnCount + j * columnCount);
        }
    }
    return array;
}

char **getTwoDArrayChar(int64_t rowCount, int64_t columnCount) {
    char **array = st_calloc(rowCount, sizeof(char *));
    array[0] = (char *) st_calloc(columnCount * rowCount, sizeof(char));
    for (int64_t i = 1; i < rowCount; i++) {
        array[i] = array[0] + i * columnCount;
    }
    return array;
}


// todo rethink this
#define MAX_TOTAL_WEIGHT 256.0

uint8_t convertTotalWeightToUInt8(double totalWeight) {
    // convert to "depth space"
    totalWeight /= PAIR_ALIGNMENT_PROB_1;
    // cap at 64
    if (totalWeight > MAX_TOTAL_WEIGHT) totalWeight = MAX_TOTAL_WEIGHT;
    // convert to uint8
    return (uint8_t) (totalWeight / MAX_TOTAL_WEIGHT * (UINT8_MAX - 1));
}

uint8_t normalizeWeightToUInt8(double totalWeight, double weight) {
    return (uint8_t) (weight / totalWeight * (UINT8_MAX - 1));
}


void
writeSimpleWeightHelenFeaturesHDF5(Alphabet *alphabet, HelenFeatureHDF5FileInfo *hdf5FileInfo, char *outputFileBase,
                                   BamChunk *bamChunk, bool outputLabels, stList *features, int64_t featureStartIdx,
                                   int64_t featureEndIdxInclusive) {

    // count features, create feature array
    uint64_t featureCount = 0;
    for (int64_t i = featureStartIdx; i <= featureEndIdxInclusive; i++) {
        PoaFeatureSimpleWeight *feature = stList_get(features, i);
        while (feature != NULL) {
            featureCount++;
            feature = feature->nextInsert;
        }
    }
    if (featureCount < HDF5_FEATURE_SIZE && outputLabels) {
        char *logIdentifier = getLogIdentifier();
        st_logInfo(" %s Feature count %"PRId64" less than minimum of %d\n", logIdentifier, featureCount,
                   HDF5_FEATURE_SIZE);
        free(logIdentifier);
        return;
    }

    // get all feature data into an array
    uint32_t **positionData = getTwoDArrayUInt32(featureCount, 2);
    int64_t columnCount = POAFEATURE_SIMPLE_WEIGHT_TOTAL_SIZE;
    uint8_t **normalizationData = getTwoDArrayUInt8(featureCount, 1);
    uint8_t **imageData = getTwoDArrayUInt8(featureCount, columnCount);
    char **labelCharacterData = NULL;
    if (outputLabels) {
        labelCharacterData = getTwoDArrayChar(featureCount, 1);
    }

    // add all data to features
    featureCount = 0;
    for (int64_t i = featureStartIdx; i <= featureEndIdxInclusive; i++) {
        PoaFeatureSimpleWeight *feature = stList_get(features, i);

        // get total weight
        double totalWeight = 0;
        for (int64_t j = 0; j < columnCount; j++) {
            totalWeight += feature->weights[j];
        }

        while (feature != NULL) {
            //positions
            positionData[featureCount][0] = (uint32_t) feature->refPosition;
            positionData[featureCount][1] = (uint32_t) feature->insertPosition;

            // normalization and image data
            normalizationData[featureCount][0] = convertTotalWeightToUInt8(totalWeight);
            for (int64_t j = 0; j < columnCount; j++) {
                imageData[featureCount][j] = normalizeWeightToUInt8(totalWeight, feature->weights[j]);
            }

            // nucleotide weights
            int64_t pos;
            for (int64_t symbol = 0; symbol < SYMBOL_NUMBER_NO_N; symbol++) {
                pos = PoaFeature_SimpleWeight_charIndex((Symbol) symbol, TRUE);
                imageData[featureCount][pos] = normalizeWeightToUInt8(totalWeight, feature->weights[pos]);
                pos = PoaFeature_SimpleWeight_charIndex((Symbol) symbol, FALSE);
                imageData[featureCount][pos] = normalizeWeightToUInt8(totalWeight, feature->weights[pos]);
            }
            // gap weights
            pos = PoaFeature_SimpleWeight_gapIndex(TRUE);
            imageData[featureCount][pos] = normalizeWeightToUInt8(totalWeight, feature->weights[pos]);
            pos = PoaFeature_SimpleWeight_gapIndex(FALSE);
            imageData[featureCount][pos] = normalizeWeightToUInt8(totalWeight, feature->weights[pos]);

            // potentially labels
            if (outputLabels) {
                Symbol label = alphabet->convertCharToSymbol(feature->label);
                labelCharacterData[featureCount][0] = (uint8_t) (alphabet->convertSymbolToChar(label) == 'N' ? 0 :
                                                                 label + 1);
            }

            // increment
            featureCount++;
            feature = feature->nextInsert;

        }
    }

    /*
     * Get hdf5 data set up
     */

    hid_t status;
    hid_t stringType = H5Tcopy(H5T_C_S1);
    status = H5Tset_size(stringType, strlen(bamChunk->refSeqName) + 1);

    // so that we can produce chunks smaller than HDF5_FEATURE_SIZE (not used during training)
    hsize_t featureSize = (hsize_t) (featureCount < HDF5_FEATURE_SIZE ? featureCount : HDF5_FEATURE_SIZE);

    hsize_t metadataDimension[1] = {1};
    hsize_t postionDimension[2] = {featureSize, 2};
    hsize_t labelCharacterDimension[2] = {featureSize, 1};
    hsize_t normalizationDimension[2] = {featureSize, 1};
    hsize_t imageDimension[2] = {featureSize, (hsize_t) columnCount};

    hid_t metadataSpace = H5Screate_simple(1, metadataDimension, NULL);
    hid_t positionSpace = H5Screate_simple(2, postionDimension, NULL);
    hid_t labelCharacterSpace = H5Screate_simple(2, labelCharacterDimension, NULL);
    hid_t normalizationSpace = H5Screate_simple(2, normalizationDimension, NULL);
    hid_t imageSpace = H5Screate_simple(2, imageDimension, NULL);

    /*
     * Write features to files
     */

    // each file must have exactly 1000 features
    int64_t totalFeatureFiles =
            (int64_t) (featureCount / HDF5_FEATURE_SIZE) + (featureCount % HDF5_FEATURE_SIZE == 0 ? 0 : 1);
    int64_t featureOffset = 0;
    if (featureCount >= HDF5_FEATURE_SIZE) {
        featureOffset = (int64_t) ((HDF5_FEATURE_SIZE * totalFeatureFiles - featureCount) /
                                   (int64_t) (featureCount / HDF5_FEATURE_SIZE));
    }

    for (int64_t featureIndex = 0; featureIndex < totalFeatureFiles; featureIndex++) {
        // get start pos
        int64_t chunkFeatureStartIdx = (HDF5_FEATURE_SIZE * featureIndex) - (featureOffset * featureIndex);
        if (featureIndex + 1 == totalFeatureFiles && featureCount >= HDF5_FEATURE_SIZE) {
            chunkFeatureStartIdx = featureCount - HDF5_FEATURE_SIZE;
        }

        // create group
        char *outputGroup = stString_print("images/%s.%"PRId64, outputFileBase, featureIndex);
        hid_t group = H5Gcreate(hdf5FileInfo->file, outputGroup, hdf5FileInfo->groupPropertyList, H5P_DEFAULT,
                                H5P_DEFAULT);

        // write metadata
        hid_t contigDataset = H5Dcreate(group, "contig", stringType, metadataSpace, H5P_DEFAULT, H5P_DEFAULT,
                                        H5P_DEFAULT);
        status |= H5Dwrite(contigDataset, stringType, H5S_ALL, H5S_ALL, H5P_DEFAULT, bamChunk->refSeqName);
        hid_t contigStartDataset = H5Dcreate(group, "contig_start", hdf5FileInfo->int64Type, metadataSpace,
                                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite(contigStartDataset, hdf5FileInfo->int64Type,
                           H5S_ALL, H5S_ALL, H5P_DEFAULT, &bamChunk->chunkBoundaryStart);
        hid_t contigEndDataset = H5Dcreate(group, "contig_end", hdf5FileInfo->int64Type, metadataSpace,
                                           H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite(contigEndDataset, hdf5FileInfo->int64Type,
                           H5S_ALL, H5S_ALL, H5P_DEFAULT, &bamChunk->chunkBoundaryEnd);
        hid_t chunkIndexDataset = H5Dcreate(group, "feature_chunk_idx", hdf5FileInfo->int64Type, metadataSpace,
                                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite(chunkIndexDataset, hdf5FileInfo->int64Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &featureIndex);

        // write position info
        hid_t positionDataset = H5Dcreate(group, "position", hdf5FileInfo->uint32Type, positionSpace,
                                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite(positionDataset, hdf5FileInfo->uint32Type,
                           H5S_ALL, H5S_ALL, H5P_DEFAULT, positionData[chunkFeatureStartIdx]);

        // write normalization data
        hid_t normalizationDataset = H5Dcreate(group, "normalization", hdf5FileInfo->uint8Type, normalizationSpace,
                                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite(normalizationDataset, hdf5FileInfo->uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                           normalizationData[chunkFeatureStartIdx]);

        // write rle data
        hid_t imageDataset = H5Dcreate(group, "image", hdf5FileInfo->uint8Type, imageSpace,
                                       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite(imageDataset, hdf5FileInfo->uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                           imageData[chunkFeatureStartIdx]);

        // if labels, add all these too
        if (outputLabels) {
            hid_t labelCharacterDataset = H5Dcreate(group, "label_base", hdf5FileInfo->uint8Type, labelCharacterSpace,
                                                    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            status |= H5Dwrite(labelCharacterDataset, hdf5FileInfo->uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                               labelCharacterData[chunkFeatureStartIdx]);
            status |= H5Dclose(labelCharacterDataset);
        }

        // cleanup
        status |= H5Dclose(contigDataset);
        status |= H5Dclose(contigStartDataset);
        status |= H5Dclose(contigEndDataset);
        status |= H5Dclose(chunkIndexDataset);
        status |= H5Dclose(positionDataset);
        status |= H5Dclose(normalizationDataset);
        status |= H5Dclose(imageDataset);
        status |= H5Gclose(group);
        free(outputGroup);
    }


    // cleanup
    free(normalizationData[0]);
    free(normalizationData);
    free(imageData[0]);
    free(imageData);
    free(positionData[0]);
    free(positionData);
    status |= H5Tclose(stringType);
    status |= H5Sclose(metadataSpace);
    status |= H5Sclose(normalizationSpace);
    status |= H5Sclose(imageSpace);
    status |= H5Sclose(positionSpace);
    status |= H5Sclose(labelCharacterSpace);
    if (outputLabels) {
        free(labelCharacterData[0]);
        free(labelCharacterData);
    }

    if (status) {
        char *logIdentifier = getLogIdentifier();
        st_logInfo(" %s Error writing HELEN features to HDF5 files: %s\n", logIdentifier, outputFileBase);
        free(logIdentifier);
    }
}

void
writeSplitRleWeightHelenFeaturesHDF5(Alphabet *alphabet, HelenFeatureHDF5FileInfo *hdf5FileInfo, char *outputFileBase,
                                     BamChunk *bamChunk,
                                     bool outputLabels, stList *features, int64_t featureStartIdx,
                                     int64_t featureEndIdxInclusive, const int64_t maxRunLength) {

    herr_t status = 0;

    /*
     * Get feature data set up
     */

    // count features, create feature array
    uint64_t featureCount = 0;
    for (int64_t i = featureStartIdx; i <= featureEndIdxInclusive; i++) {
        PoaFeatureSplitRleWeight *feature = stList_get(features, i);
        while (feature != NULL) {
            PoaFeatureSplitRleWeight *rlFeature = feature;
            while (rlFeature != NULL) {
                featureCount++;
                rlFeature = rlFeature->nextRunLength;
            }
            feature = feature->nextInsert;
        }
    }
    if (featureCount < HDF5_FEATURE_SIZE && outputLabels) {
        char *logIdentifier = getLogIdentifier();
        st_logInfo(" %s Feature count %"PRId64" less than minimum of %d\n", logIdentifier, featureCount,
                   HDF5_FEATURE_SIZE);
        free(logIdentifier);
        return;
    }

    // get all feature data into an array
    uint32_t **positionData = getTwoDArrayUInt32(featureCount, 3);
    int64_t rleNucleotideColumnCount = ((SYMBOL_NUMBER - 1) * (maxRunLength + 1) + 1) * 2;
    uint8_t **normalizationData = getTwoDArrayUInt8(featureCount, 1);
    uint8_t **imageData = getTwoDArrayUInt8(featureCount, rleNucleotideColumnCount);
    uint8_t **labelCharacterData = NULL;
    uint8_t **labelRunLengthData = NULL;
    if (outputLabels) {
        labelCharacterData = getTwoDArrayUInt8(featureCount, 1);
        labelRunLengthData = getTwoDArrayUInt8(featureCount, 1);
    }

    // add all data to features
    featureCount = 0;
    for (int64_t i = featureStartIdx; i <= featureEndIdxInclusive; i++) {
        PoaFeatureSplitRleWeight *refFeature = stList_get(features, i);

        // total weight is calculated for the very first refPos feature, used for all inserts and run lengths
        double totalWeight = 0;
        for (int64_t j = 0; j < rleNucleotideColumnCount; j++) {
            totalWeight += refFeature->weights[j];
        }

        // iterate over all insert features
        PoaFeatureSplitRleWeight *insFeature = refFeature;
        while (insFeature != NULL) {

            // iterate over all run length features
            PoaFeatureSplitRleWeight *rlFeature = insFeature;
            while (rlFeature != NULL) {
                // position
                positionData[featureCount][0] = (uint32_t) rlFeature->refPosition;
                positionData[featureCount][1] = (uint32_t) rlFeature->insertPosition;
                positionData[featureCount][2] = (uint32_t) rlFeature->runLengthPosition;

                // normalization
                normalizationData[featureCount][0] = convertTotalWeightToUInt8(totalWeight);

                // copy weights over (into normalized uint8 space)
                for (int64_t j = 0; j < rleNucleotideColumnCount; j++) {
                    imageData[featureCount][j] =
                            normalizeWeightToUInt8(totalWeight, rlFeature->weights[j]);
                }

                // labels
                if (outputLabels) {
                    Symbol label = alphabet->convertCharToSymbol(rlFeature->labelChar);
                    labelCharacterData[featureCount][0] = (uint8_t) (alphabet->convertSymbolToChar(label) == 'N' ?
                                                                     0 : label + 1);
                    labelRunLengthData[featureCount][0] = (uint8_t) (alphabet->convertSymbolToChar(label) == 'N' ?
                                                                     0 : rlFeature->labelRunLength);
                    if (labelRunLengthData[featureCount][0] > maxRunLength) {
                        st_errAbort("Encountered run length of %d (max %"PRId64") in chunk %s:%"PRId64"-%"PRId64,
                                    labelRunLengthData[featureCount][0], maxRunLength, bamChunk->refSeqName,
                                    bamChunk->chunkBoundaryStart, bamChunk->chunkBoundaryEnd);
                    }
                }

                // iterate
                featureCount++;
                rlFeature = rlFeature->nextRunLength;
            }
            insFeature = insFeature->nextInsert;
        }
    }

    /*
     * Get hdf5 data set up
     */


    // so that we can produce chunks smaller than HDF5_FEATURE_SIZE (not used during training)
    hsize_t featureSize = (hsize_t) (featureCount < HDF5_FEATURE_SIZE ? featureCount : HDF5_FEATURE_SIZE);

    hsize_t metadataDimension[1] = {1};
    hsize_t postionDimension[2] = {featureSize, 3};
    hsize_t labelCharacterDimension[2] = {featureSize, 1};
    hsize_t labelRunLengthDimension[2] = {featureSize, 1};
    hsize_t normalizationDimension[2] = {featureSize, 1};
    hsize_t imageDimension[2] = {featureSize, (hsize_t) rleNucleotideColumnCount};

    hid_t metadataSpace = H5Screate_simple(1, metadataDimension, NULL);
    hid_t positionSpace = H5Screate_simple(2, postionDimension, NULL);
    hid_t labelCharacterSpace = H5Screate_simple(2, labelCharacterDimension, NULL);
    hid_t labelRunLengthSpace = H5Screate_simple(2, labelRunLengthDimension, NULL);
    hid_t normalizationSpace = H5Screate_simple(2, normalizationDimension, NULL);
    hid_t imageSpace = H5Screate_simple(2, imageDimension, NULL);

    hid_t stringType = H5Tcopy(H5T_C_S1);
    H5Tset_size(stringType, strlen(bamChunk->refSeqName) + 1);


    /*
     * Write features to files
     */

    // each file must have exactly 1000 features
    int64_t totalFeatureFiles =
            (int64_t) (featureCount / HDF5_FEATURE_SIZE) + (featureCount % HDF5_FEATURE_SIZE == 0 ? 0 : 1);
    int64_t featureOffset = 0;
    if (featureCount >= HDF5_FEATURE_SIZE) {
        featureOffset = (int64_t) ((HDF5_FEATURE_SIZE * totalFeatureFiles - featureCount) /
                                   (int64_t) (featureCount / HDF5_FEATURE_SIZE));
    }
    for (int64_t featureIndex = 0; featureIndex < totalFeatureFiles; featureIndex++) {
        // get start pos
        int64_t chunkFeatureStartIdx = (HDF5_FEATURE_SIZE * featureIndex) - (featureOffset * featureIndex);
        if (featureIndex + 1 == totalFeatureFiles && featureCount >= HDF5_FEATURE_SIZE) {
            chunkFeatureStartIdx = featureCount - HDF5_FEATURE_SIZE;
        }

        // create group
        char *outputGroup = stString_print("images/%s.%"PRId64, outputFileBase, featureIndex);
        hid_t group = H5Gcreate(hdf5FileInfo->file, outputGroup, hdf5FileInfo->groupPropertyList, H5P_DEFAULT,
                                H5P_DEFAULT);

        // write metadata
        hid_t contigDataset = H5Dcreate(group, "contig", stringType, metadataSpace, H5P_DEFAULT, H5P_DEFAULT,
                                        H5P_DEFAULT);
        status |= H5Dwrite(contigDataset, stringType, H5S_ALL, H5S_ALL, H5P_DEFAULT, bamChunk->refSeqName);
        hid_t contigStartDataset = H5Dcreate(group, "contig_start", hdf5FileInfo->int64Type, metadataSpace, H5P_DEFAULT,
                                             H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite(contigStartDataset, hdf5FileInfo->int64Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                           &bamChunk->chunkBoundaryStart);
        hid_t contigEndDataset = H5Dcreate(group, "contig_end", hdf5FileInfo->int64Type, metadataSpace, H5P_DEFAULT,
                                           H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite(contigEndDataset, hdf5FileInfo->int64Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                           &bamChunk->chunkBoundaryEnd);
        hid_t chunkIndexDataset = H5Dcreate(group, "feature_chunk_idx", hdf5FileInfo->int64Type, metadataSpace,
                                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite(chunkIndexDataset, hdf5FileInfo->int64Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &featureIndex);

        // write position info
        hid_t positionDataset = H5Dcreate(group, "position", hdf5FileInfo->uint32Type, positionSpace, H5P_DEFAULT,
                                          H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite(positionDataset, hdf5FileInfo->uint32Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                           positionData[chunkFeatureStartIdx]);

        // write rle data
        hid_t imageDataset = H5Dcreate(group, "image", hdf5FileInfo->uint8Type, imageSpace, H5P_DEFAULT, H5P_DEFAULT,
                                       H5P_DEFAULT);
        status |= H5Dwrite(imageDataset, hdf5FileInfo->uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                           imageData[chunkFeatureStartIdx]);
        hid_t normalizationDataset = H5Dcreate(group, "normalization", hdf5FileInfo->uint8Type, normalizationSpace,
                                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite(normalizationDataset, hdf5FileInfo->uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                           normalizationData[chunkFeatureStartIdx]);

        // if labels, add all these too
        if (outputLabels) {
            hid_t labelCharacterDataset = H5Dcreate(group, "label_base", hdf5FileInfo->uint8Type, labelCharacterSpace,
                                                    H5P_DEFAULT,
                                                    H5P_DEFAULT, H5P_DEFAULT);
            status |= H5Dwrite(labelCharacterDataset, hdf5FileInfo->uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                               labelCharacterData[chunkFeatureStartIdx]);
            hid_t labelRunLengthDataset = H5Dcreate(group, "label_run_length", hdf5FileInfo->uint8Type,
                                                    labelRunLengthSpace,
                                                    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            status |= H5Dwrite(labelRunLengthDataset, hdf5FileInfo->uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                               labelRunLengthData[chunkFeatureStartIdx]);

            status |= H5Dclose(labelCharacterDataset);
            status |= H5Dclose(labelRunLengthDataset);
        }

        // cleanup
        status |= H5Dclose(contigDataset);
        status |= H5Dclose(contigStartDataset);
        status |= H5Dclose(contigEndDataset);
        status |= H5Dclose(chunkIndexDataset);
        status |= H5Dclose(positionDataset);
        status |= H5Dclose(imageDataset);
        status |= H5Dclose(normalizationDataset);
        status |= H5Gclose(group);
        free(outputGroup);
    }

    // cleanup
    free(imageData[0]);
    free(imageData);
    free(normalizationData[0]);
    free(normalizationData);
    free(positionData[0]);
    free(positionData);
    status |= H5Sclose(metadataSpace);
    status |= H5Sclose(positionSpace);
    status |= H5Sclose(imageSpace);
    status |= H5Sclose(normalizationSpace);
    status |= H5Sclose(labelRunLengthSpace);
    status |= H5Sclose(labelCharacterSpace);
    status |= H5Tclose(stringType);
    if (outputLabels) {
        free(labelCharacterData[0]);
        free(labelCharacterData);
        free(labelRunLengthData[0]);
        free(labelRunLengthData);
    }

    if (status) {
        char *logIdentifier = getLogIdentifier();
        st_logInfo(" %s Error writing HELEN features to HDF5 files: %s\n", logIdentifier, outputFileBase);
        free(logIdentifier);
    }
}


void
writeChannelRleWeightHelenFeaturesHDF5(Alphabet *alphabet, HelenFeatureHDF5FileInfo *hdf5FileInfo, char *outputFileBase,
                                       BamChunk *bamChunk, bool outputLabels, stList *features, int64_t featureStartIdx,
                                       int64_t featureEndIdxInclusive, const int64_t maxRunLength) {

    herr_t status = 0;

    /*
     * Get feature data set up
     */

    // count features, create feature array
    uint64_t featureCount = 0;
    for (int64_t i = featureStartIdx; i <= featureEndIdxInclusive; i++) {
        PoaFeatureChannelRleWeight *feature = stList_get(features, i);
        while (feature != NULL) {
            PoaFeatureChannelRleWeight *rlFeature = feature;
            while (rlFeature != NULL) {
                featureCount++;
                rlFeature = rlFeature->nextRunLength;
            }
            feature = feature->nextInsert;
        }
    }
    if (featureCount < HDF5_FEATURE_SIZE && outputLabels) {
        char *logIdentifier = getLogIdentifier();
        st_logInfo(" %s Feature count %"PRId64" less than minimum of %d\n", logIdentifier, featureCount,
                   HDF5_FEATURE_SIZE);
        free(logIdentifier);
        return;
    }

    // sizes
    int64_t nucleotideColumnCount = SYMBOL_NUMBER * 2; //ACGTGap x {fwd,bwd}
    int64_t runLengthColumnCount = (maxRunLength + 1) * 2; //(runLenght + 0) x {fwd,bwd}

    // get all feature data into an array
    uint32_t **positionData = getTwoDArrayUInt32(featureCount, 3);
    uint8_t **normalizationData = getTwoDArrayUInt8(featureCount, 1);
    uint8_t **nucleotideData = getTwoDArrayUInt8(featureCount, nucleotideColumnCount);
    uint8_t ***runLengthData = getThreeDArrayUInt8(featureCount, runLengthColumnCount, (SYMBOL_NUMBER - 1));
    //uint8_t ***runLengthData = getThreeDArrayUInt8(featureCount, (SYMBOL_NUMBER - 1), runLengthColumnCount);
    uint8_t **labelCharacterData = NULL;
    uint8_t **labelRunLengthData = NULL;
    if (outputLabels) {
        labelCharacterData = getTwoDArrayUInt8(featureCount, 1);
        labelRunLengthData = getTwoDArrayUInt8(featureCount, 1);
    }

    // add all data to features
    featureCount = 0;
    for (int64_t i = featureStartIdx; i <= featureEndIdxInclusive; i++) {
        PoaFeatureChannelRleWeight *refFeature = stList_get(features, i);

        // total weight is calculated for the very first refPos feature, used for all inserts and run lengths
        double totalWeight = 0;
        for (int64_t j = 0; j < nucleotideColumnCount; j++) {
            totalWeight += refFeature->nucleotideWeights[j];
        }

        // iterate over all insert features
        PoaFeatureChannelRleWeight *insFeature = refFeature;
        while (insFeature != NULL) {

            // iterate over all run length features
            PoaFeatureChannelRleWeight *rlFeature = insFeature;
            while (rlFeature != NULL) {
                // position
                positionData[featureCount][0] = (uint32_t) rlFeature->refPosition;
                positionData[featureCount][1] = (uint32_t) rlFeature->insertPosition;
                positionData[featureCount][2] = (uint32_t) rlFeature->runLengthPosition;

                // normalization
                normalizationData[featureCount][0] = convertTotalWeightToUInt8(totalWeight);

                // copy weights over (into normalized uint8 space)
                for (int64_t c = 0; c < SYMBOL_NUMBER - 1; c++) {
                    // overall nucl count
                    nucleotideData[featureCount][c * 2 + POS_STRAND_IDX] = normalizeWeightToUInt8(totalWeight,
                                                                                                  rlFeature->nucleotideWeights[PoaFeature_ChannelRleWeight_charNuclIndex(
                                                                                                          c, TRUE)]);
                    nucleotideData[featureCount][c * 2 + NEG_STRAND_IDX] = normalizeWeightToUInt8(totalWeight,
                                                                                                  rlFeature->nucleotideWeights[PoaFeature_ChannelRleWeight_charNuclIndex(
                                                                                                          c, FALSE)]);

                    // run length counts
                    for (int64_t r = 0; r <= maxRunLength; r++) {
                        runLengthData[featureCount][r * 2 + POS_STRAND_IDX][c] =
                                normalizeWeightToUInt8(totalWeight, rlFeature->runLengthWeights[
                                        PoaFeature_ChannelRleWeight_charRLIndex(maxRunLength, c, r, TRUE)]);
                        runLengthData[featureCount][r * 2 + NEG_STRAND_IDX][c] =
                                normalizeWeightToUInt8(totalWeight, rlFeature->runLengthWeights[
                                        PoaFeature_ChannelRleWeight_charRLIndex(maxRunLength, c, r, FALSE)]);

                        //runLengthData[featureCount][c][r * 2 + POS_STRAND_IDX] =
                        //        normalizeWeightToUInt8(totalWeight, rlFeature->runLengthWeights[
                        //                PoaFeature_ChannelRleWeight_charRLIndex(maxRunLength, c, r, TRUE)]);
                        //runLengthData[featureCount][c][r * 2 + NEG_STRAND_IDX] =
                        //        normalizeWeightToUInt8(totalWeight, rlFeature->runLengthWeights[
                        //                PoaFeature_ChannelRleWeight_charRLIndex(maxRunLength, c, r, FALSE)]);
                    }
                }
                // gap counts
                nucleotideData[featureCount][SYMBOL_NUMBER_NO_N * 2 + 0 + POS_STRAND_IDX] = normalizeWeightToUInt8(
                        totalWeight, rlFeature->nucleotideWeights[PoaFeature_ChannelRleWeight_gapNuclIndex(TRUE)]);
                nucleotideData[featureCount][SYMBOL_NUMBER_NO_N * 2 + NEG_STRAND_IDX] = normalizeWeightToUInt8(
                        totalWeight, rlFeature->nucleotideWeights[PoaFeature_ChannelRleWeight_gapNuclIndex(FALSE)]);

                // labels
                if (outputLabels) {
                    Symbol label = alphabet->convertCharToSymbol(rlFeature->labelChar);
                    labelCharacterData[featureCount][0] = (uint8_t) (alphabet->convertSymbolToChar(label) == 'N' ?
                                                                     0 : label + 1);
                    labelRunLengthData[featureCount][0] = (uint8_t) (alphabet->convertSymbolToChar(label) == 'N' ?
                                                                     0 : rlFeature->labelRunLength);
                    if (labelRunLengthData[featureCount][0] > maxRunLength) {
                        st_errAbort("Encountered run length of %d (max %"PRId64") in chunk %s:%"PRId64"-%"PRId64,
                                    labelRunLengthData[featureCount][0], maxRunLength, bamChunk->refSeqName,
                                    bamChunk->chunkBoundaryStart, bamChunk->chunkBoundaryEnd);
                    }
                }

                // iterate
                featureCount++;
                rlFeature = rlFeature->nextRunLength;
            }
            insFeature = insFeature->nextInsert;
        }
    }

    /*
     * Get hdf5 data set up
     */


    // so that we can produce chunks smaller than HDF5_FEATURE_SIZE (not used during training)
    hsize_t featureSize = (hsize_t) (featureCount < HDF5_FEATURE_SIZE ? featureCount : HDF5_FEATURE_SIZE);

    hsize_t metadataDimension[1] = {1};
    hsize_t postionDimension[2] = {featureSize, 3};
    hsize_t labelCharacterDimension[2] = {featureSize, 1};
    hsize_t labelRunLengthDimension[2] = {featureSize, 1};
    hsize_t normalizationDimension[2] = {featureSize, 1};
    hsize_t nucleotideDimension[2] = {featureSize, (hsize_t) nucleotideColumnCount};
    hsize_t runLengthDimension[3] = {featureSize, (hsize_t) runLengthColumnCount, (hsize_t) SYMBOL_NUMBER - 1};
//    hsize_t runLengthDimension[3] = {featureSize, (hsize_t) SYMBOL_NUMBER - 1, (hsize_t) runLengthColumnCount};

    hid_t metadataSpace = H5Screate_simple(1, metadataDimension, NULL);
    hid_t positionSpace = H5Screate_simple(2, postionDimension, NULL);
    hid_t labelCharacterSpace = H5Screate_simple(2, labelCharacterDimension, NULL);
    hid_t labelRunLengthSpace = H5Screate_simple(2, labelRunLengthDimension, NULL);
    hid_t normalizationSpace = H5Screate_simple(2, normalizationDimension, NULL);
    hid_t nucleotideSpace = H5Screate_simple(2, nucleotideDimension, NULL);
    hid_t runLengthSpace = H5Screate_simple(3, runLengthDimension, NULL);

    hid_t stringType = H5Tcopy(H5T_C_S1);
    H5Tset_size(stringType, strlen(bamChunk->refSeqName) + 1);


    /*
     * Write features to files
     */

    // each file must have exactly 1000 features
    int64_t totalFeatureFiles =
            (int64_t) (featureCount / HDF5_FEATURE_SIZE) + (featureCount % HDF5_FEATURE_SIZE == 0 ? 0 : 1);
    int64_t featureOffset = 0;
    if (featureCount >= HDF5_FEATURE_SIZE) {
        featureOffset = (int64_t) ((HDF5_FEATURE_SIZE * totalFeatureFiles - featureCount) /
                                   (int64_t) (featureCount / HDF5_FEATURE_SIZE));
    }
    for (int64_t featureIndex = 0; featureIndex < totalFeatureFiles; featureIndex++) {
        // get start pos
        int64_t chunkFeatureStartIdx = (HDF5_FEATURE_SIZE * featureIndex) - (featureOffset * featureIndex);
        if (featureIndex + 1 == totalFeatureFiles && featureCount >= HDF5_FEATURE_SIZE) {
            chunkFeatureStartIdx = featureCount - HDF5_FEATURE_SIZE;
        }

        // create group
        char *outputGroup = stString_print("images/%s.%"PRId64, outputFileBase, featureIndex);
        hid_t group = H5Gcreate(hdf5FileInfo->file, outputGroup, hdf5FileInfo->groupPropertyList, H5P_DEFAULT,
                                H5P_DEFAULT);

        // write metadata
        hid_t contigDataset = H5Dcreate(group, "contig", stringType, metadataSpace, H5P_DEFAULT, H5P_DEFAULT,
                                        H5P_DEFAULT);
        status |= H5Dwrite(contigDataset, stringType, H5S_ALL, H5S_ALL, H5P_DEFAULT, bamChunk->refSeqName);
        hid_t contigStartDataset = H5Dcreate(group, "contig_start", hdf5FileInfo->int64Type, metadataSpace, H5P_DEFAULT,
                                             H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite(contigStartDataset, hdf5FileInfo->int64Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                           &bamChunk->chunkBoundaryStart);
        hid_t contigEndDataset = H5Dcreate(group, "contig_end", hdf5FileInfo->int64Type, metadataSpace, H5P_DEFAULT,
                                           H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite(contigEndDataset, hdf5FileInfo->int64Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                           &bamChunk->chunkBoundaryEnd);
        hid_t chunkIndexDataset = H5Dcreate(group, "feature_chunk_idx", hdf5FileInfo->int64Type, metadataSpace,
                                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite(chunkIndexDataset, hdf5FileInfo->int64Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &featureIndex);

        // write position info
        hid_t positionDataset = H5Dcreate(group, "position", hdf5FileInfo->uint32Type, positionSpace, H5P_DEFAULT,
                                          H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite(positionDataset, hdf5FileInfo->uint32Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                           positionData[chunkFeatureStartIdx]);

        // write nucl, rl, and norm data
        hid_t nucleotideDataset = H5Dcreate(group, "nucleotide", hdf5FileInfo->uint8Type, nucleotideSpace, H5P_DEFAULT,
                                            H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite(nucleotideDataset, hdf5FileInfo->uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                           nucleotideData[chunkFeatureStartIdx]);
        hid_t runLengthDataset = H5Dcreate(group, "runLengths", hdf5FileInfo->uint8Type, runLengthSpace, H5P_DEFAULT,
                                           H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite(runLengthDataset, hdf5FileInfo->uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                           runLengthData[chunkFeatureStartIdx][0]);
        hid_t normalizationDataset = H5Dcreate(group, "normalization", hdf5FileInfo->uint8Type, normalizationSpace,
                                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite(normalizationDataset, hdf5FileInfo->uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                           normalizationData[chunkFeatureStartIdx]);

        // if labels, add all these too
        if (outputLabels) {
            hid_t labelCharacterDataset = H5Dcreate(group, "label_base", hdf5FileInfo->uint8Type, labelCharacterSpace,
                                                    H5P_DEFAULT,
                                                    H5P_DEFAULT, H5P_DEFAULT);
            status |= H5Dwrite(labelCharacterDataset, hdf5FileInfo->uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                               labelCharacterData[chunkFeatureStartIdx]);
            hid_t labelRunLengthDataset = H5Dcreate(group, "label_run_length", hdf5FileInfo->uint8Type,
                                                    labelRunLengthSpace,
                                                    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            status |= H5Dwrite(labelRunLengthDataset, hdf5FileInfo->uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                               labelRunLengthData[chunkFeatureStartIdx]);

            status |= H5Dclose(labelCharacterDataset);
            status |= H5Dclose(labelRunLengthDataset);
        }

        // cleanup
        status |= H5Dclose(contigDataset);
        status |= H5Dclose(contigStartDataset);
        status |= H5Dclose(contigEndDataset);
        status |= H5Dclose(chunkIndexDataset);
        status |= H5Dclose(positionDataset);
        status |= H5Dclose(nucleotideDataset);
        status |= H5Dclose(runLengthDataset);
        status |= H5Dclose(normalizationDataset);
        status |= H5Gclose(group);
        free(outputGroup);
    }

    // cleanup
    free(nucleotideData[0]);
    free(nucleotideData);
    free(runLengthData[0][0]);
    free(runLengthData[0]);
    free(runLengthData);
    free(normalizationData[0]);
    free(normalizationData);
    free(positionData[0]);
    free(positionData);
    status |= H5Sclose(metadataSpace);
    status |= H5Sclose(positionSpace);
    status |= H5Sclose(nucleotideSpace);
    status |= H5Sclose(runLengthSpace);
    status |= H5Sclose(normalizationSpace);
    status |= H5Sclose(labelRunLengthSpace);
    status |= H5Sclose(labelCharacterSpace);
    status |= H5Tclose(stringType);
    if (outputLabels) {
        free(labelCharacterData[0]);
        free(labelCharacterData);
        free(labelRunLengthData[0]);
        free(labelRunLengthData);
    }

    if (status) {
        char *logIdentifier = getLogIdentifier();
        st_logInfo(" %s Error writing HELEN features to HDF5 files: %s\n", logIdentifier, outputFileBase);
        free(logIdentifier);
    }
}

void writeDiploidRleWeightHelenFeaturesHDF5(Alphabet *alphabet, HelenFeatureHDF5FileInfo* hdf5FileInfo,
        char *outputFileBase, BamChunk *bamChunk, int64_t haplotype, int64_t truthIdx,
        stList *features, int64_t featureStartIdx, int64_t featureEndIdxInclusive,
        const int64_t maxRunLength) {

    bool outputLabels = truthIdx >= 0;
    if (!outputLabels) truthIdx = 0;
    herr_t status = 0;
    /*
     * Get feature data set up
     */

    // count features, create feature array
    uint64_t featureCount = 0;
    for (int64_t i = featureStartIdx; i <= featureEndIdxInclusive; i++) {
        PoaFeatureDiploidRleWeight *feature = stList_get(features, i);
        while (feature != NULL) {
            PoaFeatureDiploidRleWeight *rlFeature = feature;
            while (rlFeature != NULL) {
                featureCount++;
                rlFeature = rlFeature->nextRunLength;
            }
            feature = feature->nextInsert;
        }
    }

    // don't write small feature sets when training
    if (featureCount < HDF5_FEATURE_SIZE && outputLabels) {
        char *logIdentifier = getLogIdentifier();
        st_logInfo(" %s Feature count %"PRId64" (H%"PRId64", T%"PRId64") less than minimum of %d\n", logIdentifier,
                featureCount, haplotype, truthIdx, HDF5_FEATURE_SIZE);
        free(logIdentifier);
        return;
    }

    // get all feature data into an array
    uint32_t **positionData = getTwoDArrayUInt32(featureCount, 3);
    int64_t rleNucleotideStrandColumnCount_singleHap = ((SYMBOL_NUMBER - 1) * (maxRunLength + 1) + 1) * 2; //(nucl*(rl+0)+gap)*strand
    uint8_t **normalizationData = getTwoDArrayUInt8(featureCount, 1);
    uint8_t **imageData = getTwoDArrayUInt8(featureCount, rleNucleotideStrandColumnCount_singleHap * 2);
    uint8_t **labelCharacterData = NULL;
    uint8_t **labelRunLengthData = NULL;
    if (outputLabels) {
        labelCharacterData = getTwoDArrayUInt8(featureCount, 1);
        labelRunLengthData = getTwoDArrayUInt8(featureCount, 1);
    }

    /*
     *   Set up feature data in HDF5 format
     */

    // add feature data H1
    featureCount = 0;
    for (int64_t i = featureStartIdx; i <= featureEndIdxInclusive; i++) {
        PoaFeatureDiploidRleWeight *refFeature = stList_get(features, i);

        // total weight is calculated for the very first refPos feature, used for all inserts and run lengths
        double totalWeightOn = 0;
        double totalWeightOff = 0;
        for (int64_t j = 0; j < rleNucleotideStrandColumnCount_singleHap; j++) {
            totalWeightOn += refFeature->weightsHOn[j];
            totalWeightOff += refFeature->weightsHOff[j];
        }

        // iterate over all insert features
        PoaFeatureDiploidRleWeight *insFeature = refFeature;
        while (insFeature != NULL) {

            // iterate over all run length features
            PoaFeatureDiploidRleWeight *rlFeature = insFeature;
            while (rlFeature != NULL) {
                // position
                positionData[featureCount][0] = (uint32_t) rlFeature->refPosition;
                positionData[featureCount][1] = (uint32_t) rlFeature->insertPosition;
                positionData[featureCount][2] = (uint32_t) rlFeature->runLengthPosition;

                // normalization
                normalizationData[featureCount][0] = convertTotalWeightToUInt8(totalWeightOn + totalWeightOff);

                // copy weights over (into normalized uint8 space)
                for (int64_t j = 0; j < rleNucleotideStrandColumnCount_singleHap; j++) {
                    imageData[featureCount][j] = normalizeWeightToUInt8(totalWeightOn, rlFeature->weightsHOn[j]);
                    imageData[featureCount][j+rleNucleotideStrandColumnCount_singleHap] =
                            normalizeWeightToUInt8(totalWeightOff, rlFeature->weightsHOff[j]);
                }

                // labels
                if (outputLabels) {
                    Symbol label = alphabet->convertCharToSymbol(rlFeature->labelChar);
                    labelCharacterData[featureCount][0] = (uint8_t) (alphabet->convertSymbolToChar(label) == 'N' ?
                                                                     0 : label + 1);
                    labelRunLengthData[featureCount][0] = (uint8_t) (alphabet->convertSymbolToChar(label) == 'N' ?
                                                                     0 : rlFeature->labelRunLength);
                    if (labelRunLengthData[featureCount][0] > maxRunLength) {
                        st_errAbort("Encountered run length of %d (max %"PRId64") in H1 chunk %s:%"PRId64"-%"PRId64,
                                    labelRunLengthData[featureCount][0], maxRunLength, bamChunk->refSeqName,
                                    bamChunk->chunkBoundaryStart, bamChunk->chunkBoundaryEnd);
                    }
                }

                // iterate
                featureCount++;
                rlFeature = rlFeature->nextRunLength;
            }
            insFeature = insFeature->nextInsert;
        }
    }

    /*
     * Get hdf5 data set up
     */
    // so that we can produce chunks smaller than HDF5_FEATURE_SIZE (not used during training)
    hsize_t featureSize = (hsize_t) (featureCount < HDF5_FEATURE_SIZE ? featureCount : HDF5_FEATURE_SIZE);

    hsize_t metadataDimension[1] = {1};
    hsize_t postionDimension[2] = {featureSize, 3};
    hsize_t labelCharacterDimension[2] = {featureSize, 1};
    hsize_t labelRunLengthDimension[2] = {featureSize, 1};
    hsize_t normalizationDimension[2] = {featureSize, 1};
    hsize_t imageDimension[2] = {featureSize, (hsize_t) rleNucleotideStrandColumnCount_singleHap * 2};

    hid_t metadataSpace = H5Screate_simple(1, metadataDimension, NULL);
    hid_t positionSpace = H5Screate_simple(2, postionDimension, NULL);
    hid_t labelCharacterSpace = H5Screate_simple(2, labelCharacterDimension, NULL);
    hid_t labelRunLengthSpace = H5Screate_simple(2, labelRunLengthDimension, NULL);
    hid_t normalizationSpace = H5Screate_simple(2, normalizationDimension, NULL);
    hid_t imageSpace = H5Screate_simple(2, imageDimension, NULL);

    hid_t stringType = H5Tcopy(H5T_C_S1);
    H5Tset_size(stringType, strlen(bamChunk->refSeqName) + 1);

    /*
     *   Write features to files
     */

    // Write for hap1
    int64_t totalFeatureFilesH1 =
            (int64_t) (featureCount / HDF5_FEATURE_SIZE) + (featureCount % HDF5_FEATURE_SIZE == 0 ? 0 : 1);
    int64_t featureOffsetH1 = 0;
    if (featureCount >= HDF5_FEATURE_SIZE) {
        featureOffsetH1 = (int64_t) ((HDF5_FEATURE_SIZE * totalFeatureFilesH1 - featureCount) /
                                   (int64_t) (featureCount / HDF5_FEATURE_SIZE));
    }
    for (int64_t featureIndex = 0; featureIndex < totalFeatureFilesH1; featureIndex++) {
        // get start pos
        int64_t chunkFeatureStartIdx = (HDF5_FEATURE_SIZE * featureIndex) - (featureOffsetH1 * featureIndex);
        if (featureIndex + 1 == totalFeatureFilesH1 && featureCount >= HDF5_FEATURE_SIZE) {
            chunkFeatureStartIdx = featureCount - HDF5_FEATURE_SIZE;
        }

        // create group
        char *outputGroup = stString_print("images/%s.H%"PRId64".U%"PRId64".F%"PRId64, outputFileBase, haplotype,
                truthIdx < 0 ? 0 : truthIdx, featureIndex);
        hid_t group = H5Gcreate(hdf5FileInfo->file, outputGroup, hdf5FileInfo->groupPropertyList, H5P_DEFAULT,
                                H5P_DEFAULT);

        // write metadata
        hid_t contigDataset = H5Dcreate(group, "contig", stringType, metadataSpace, H5P_DEFAULT, H5P_DEFAULT,
                                        H5P_DEFAULT);
        status |= H5Dwrite(contigDataset, stringType, H5S_ALL, H5S_ALL, H5P_DEFAULT, bamChunk->refSeqName);
        hid_t contigStartDataset = H5Dcreate(group, "contig_start", hdf5FileInfo->int64Type, metadataSpace, H5P_DEFAULT,
                                             H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite(contigStartDataset, hdf5FileInfo->int64Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                           &bamChunk->chunkBoundaryStart);
        hid_t contigEndDataset = H5Dcreate(group, "contig_end", hdf5FileInfo->int64Type, metadataSpace, H5P_DEFAULT,
                                           H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite(contigEndDataset, hdf5FileInfo->int64Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                           &bamChunk->chunkBoundaryEnd);
        hid_t chunkIndexDataset = H5Dcreate(group, "feature_chunk_idx", hdf5FileInfo->int64Type, metadataSpace,
                                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite(chunkIndexDataset, hdf5FileInfo->int64Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &featureIndex);
        hid_t haplotypeIndexDataset = H5Dcreate(group, "haplotype", hdf5FileInfo->int64Type, metadataSpace,
                                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite(haplotypeIndexDataset, hdf5FileInfo->int64Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &haplotype);
        hid_t truthIndexDataset = H5Dcreate(group, "truth_idx", hdf5FileInfo->int64Type, metadataSpace,
                                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite(truthIndexDataset, hdf5FileInfo->int64Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &truthIdx);

        // write position info
        hid_t positionDataset = H5Dcreate(group, "position", hdf5FileInfo->uint32Type, positionSpace, H5P_DEFAULT,
                                          H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite(positionDataset, hdf5FileInfo->uint32Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                           positionData[chunkFeatureStartIdx]);

        // write rle data
        hid_t imageDataset = H5Dcreate(group, "image", hdf5FileInfo->uint8Type, imageSpace, H5P_DEFAULT, H5P_DEFAULT,
                                       H5P_DEFAULT);
        status |= H5Dwrite(imageDataset, hdf5FileInfo->uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                           imageData[chunkFeatureStartIdx]);
        hid_t normalizationDataset = H5Dcreate(group, "normalization", hdf5FileInfo->uint8Type, normalizationSpace,
                                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite(normalizationDataset, hdf5FileInfo->uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                           normalizationData[chunkFeatureStartIdx]);

        // if labels, add all these too
        if (outputLabels) {
            hid_t labelCharacterDataset = H5Dcreate(group, "label_base", hdf5FileInfo->uint8Type, labelCharacterSpace,
                                                    H5P_DEFAULT,
                                                    H5P_DEFAULT, H5P_DEFAULT);
            status |= H5Dwrite(labelCharacterDataset, hdf5FileInfo->uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                               labelCharacterData[chunkFeatureStartIdx]);
            hid_t labelRunLengthDataset = H5Dcreate(group, "label_run_length", hdf5FileInfo->uint8Type,
                                                    labelRunLengthSpace,
                                                    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            status |= H5Dwrite(labelRunLengthDataset, hdf5FileInfo->uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                               labelRunLengthData[chunkFeatureStartIdx]);

            status |= H5Dclose(labelCharacterDataset);
            status |= H5Dclose(labelRunLengthDataset);
        }

        // cleanup
        status |= H5Dclose(contigDataset);
        status |= H5Dclose(contigStartDataset);
        status |= H5Dclose(contigEndDataset);
        status |= H5Dclose(chunkIndexDataset);
        status |= H5Dclose(haplotypeIndexDataset);
        status |= H5Dclose(truthIndexDataset);
        status |= H5Dclose(positionDataset);
        status |= H5Dclose(imageDataset);
        status |= H5Dclose(normalizationDataset);
        status |= H5Gclose(group);
        free(outputGroup);
    }

    // cleanup
    free(imageData[0]);
    free(imageData);
    free(normalizationData[0]);
    free(normalizationData);
    free(positionData[0]);
    free(positionData);
    status |= H5Sclose(metadataSpace);
    status |= H5Sclose(positionSpace);
    status |= H5Sclose(imageSpace);
    status |= H5Sclose(normalizationSpace);
    status |= H5Sclose(labelRunLengthSpace);
    status |= H5Sclose(labelCharacterSpace);
    status |= H5Tclose(stringType);
    if (outputLabels) {
        free(labelCharacterData[0]);
        free(labelCharacterData);
        free(labelRunLengthData[0]);
        free(labelRunLengthData);
    }

    if (status) {
        char *logIdentifier = getLogIdentifier();
        st_logInfo(" %s Error writing HELEN features to HDF5 files: %s\n", logIdentifier, outputFileBase);
        free(logIdentifier);
    }
}


HelenFeatureHDF5FileInfo *HelenFeatureHDF5FileInfo_construct(char *filename) {
    HelenFeatureHDF5FileInfo *fileInfo = st_calloc(1, sizeof(HelenFeatureHDF5FileInfo));
    fileInfo->filename = stString_copy(filename);
    fileInfo->file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    fileInfo->int64Type = H5Tcopy(H5T_NATIVE_UINT32);
    H5Tset_order(fileInfo->int64Type, H5T_ORDER_LE);
    fileInfo->uint32Type = H5Tcopy(H5T_NATIVE_UINT32);
    H5Tset_order(fileInfo->uint32Type, H5T_ORDER_LE);
    fileInfo->uint8Type = H5Tcopy(H5T_NATIVE_UINT8);
    H5Tset_order(fileInfo->uint8Type, H5T_ORDER_LE);
    fileInfo->floatType = H5Tcopy(H5T_NATIVE_FLOAT);
    H5Tset_order(fileInfo->floatType, H5T_ORDER_LE);
    fileInfo->groupPropertyList = H5Pcreate(H5P_LINK_CREATE);
    H5Pset_create_intermediate_group(fileInfo->groupPropertyList, 1);
    return fileInfo;
}

void HelenFeatureHDF5FileInfo_destruct(HelenFeatureHDF5FileInfo *fileInfo) {
    free(fileInfo->filename);
    H5Tclose(fileInfo->int64Type);
    H5Tclose(fileInfo->uint32Type);
    H5Tclose(fileInfo->uint8Type);
    H5Tclose(fileInfo->floatType);
    H5Pclose(fileInfo->groupPropertyList);
    H5Fclose(fileInfo->file);
    free(fileInfo);
}

HelenFeatureHDF5FileInfo **openHelenFeatureHDF5FilesByThreadCount(char *filenameBase, int64_t threadCount) {
    HelenFeatureHDF5FileInfo **infoArray = st_calloc(threadCount, sizeof(HelenFeatureHDF5FileInfo *));
    for (int64_t i = 0; i < threadCount; i++) {
        char *filename = stString_print("%s.T%02"PRId64".h5", filenameBase, i);
        infoArray[i] = HelenFeatureHDF5FileInfo_construct(filename);
        free(filename);
    }
    return infoArray;
}

#endif