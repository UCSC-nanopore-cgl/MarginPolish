/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com) & Arthur Rand (arand@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <sys/stat.h>
#include "margin.h"

/*
 * Get model parameters from params file.
 * Set hmm parameters.
*/

stRPHmmParameters *stRPHmmParameters_construct() {
    // Params object
    stRPHmmParameters *params = st_calloc(1, sizeof(stRPHmmParameters));

    // More variables for hmm stuff
    params->maxCoverageDepth = MAX_READ_PARTITIONING_DEPTH;
    params->maxNotSumTransitions = true;
    params->minPartitionsInAColumn = 50;
    params->maxPartitionsInAColumn = 200;
    params->minPosteriorProbabilityForPartition = 0.001;
    params->minReadCoverageToSupportPhasingBetweenHeterozygousSites = 0;

    // Other marginPhase program options
    params->roundsOfIterativeRefinement = 0;
    params->includeInvertedPartitions = true;
    params->includeAncestorSubProb = true;
    params->minPhredScoreForHaplotypePartition = 0;
    params->stitchWithPrimaryReadsOnly = TRUE;

    // bubble parameters
    params->includeHomozygousVCFEntries = FALSE;
    params->onlyUsePassVCFEntries = TRUE;
    params->onlyUseSNPVCFEntries = FALSE;
    params->useVariantSelectionAdaptiveSampling = TRUE;
    params->variantSelectionAdaptiveSamplingPrimaryThreshold = 0.9;
    params->variantSelectionAdaptiveSamplingDesiredBasepairsPerVariant = 1000;
    params->minVariantQuality = 0.15;
    params->updateAllOutputVCFFormatFields = TRUE;
    params->phasesetMinBinomialReadSplitLikelihood = .0001;
    params->phasesetMaxDiscordantRatio = .1;
    params->bubbleFindingIterations = 1;
    params->bubbleMinBinomialStrandLikelihood = .05;
    params->bubbleMinBinomialReadSplitLikelihood = .05;


    return params;
}

stRPHmmParameters *stRPHmmParameters_copy(stRPHmmParameters *toCopy) {
    // Params object
    stRPHmmParameters *params = st_calloc(1, sizeof(stRPHmmParameters));

    // More variables for hmm stuff
    params->maxCoverageDepth = toCopy->maxCoverageDepth;
    params->maxNotSumTransitions = toCopy->maxNotSumTransitions;
    params->minPartitionsInAColumn = toCopy->minPartitionsInAColumn;
    params->maxPartitionsInAColumn = toCopy->maxPartitionsInAColumn;
    params->minPosteriorProbabilityForPartition = toCopy->minPosteriorProbabilityForPartition;
    params->minReadCoverageToSupportPhasingBetweenHeterozygousSites = toCopy->minReadCoverageToSupportPhasingBetweenHeterozygousSites;

    // Other marginPhase program options
    params->roundsOfIterativeRefinement = toCopy->roundsOfIterativeRefinement;
    params->includeInvertedPartitions = toCopy->includeInvertedPartitions;
    params->includeAncestorSubProb = toCopy->includeAncestorSubProb;
    params->minPhredScoreForHaplotypePartition = toCopy->minPhredScoreForHaplotypePartition;
    params->stitchWithPrimaryReadsOnly = toCopy->stitchWithPrimaryReadsOnly;

    // bubble parameters
    params->includeHomozygousVCFEntries = toCopy->includeHomozygousVCFEntries;
    params->onlyUsePassVCFEntries = toCopy->onlyUsePassVCFEntries;
    params->onlyUseSNPVCFEntries = toCopy->onlyUseSNPVCFEntries;
    params->useVariantSelectionAdaptiveSampling = toCopy->useVariantSelectionAdaptiveSampling;
    params->variantSelectionAdaptiveSamplingPrimaryThreshold = toCopy->variantSelectionAdaptiveSamplingPrimaryThreshold;
    params->variantSelectionAdaptiveSamplingDesiredBasepairsPerVariant = toCopy->variantSelectionAdaptiveSamplingDesiredBasepairsPerVariant;
    params->minVariantQuality = toCopy->minVariantQuality;
    params->updateAllOutputVCFFormatFields = toCopy->updateAllOutputVCFFormatFields;
    params->phasesetMinBinomialReadSplitLikelihood = toCopy->phasesetMinBinomialReadSplitLikelihood;
    params->phasesetMaxDiscordantRatio = toCopy->phasesetMaxDiscordantRatio;
    params->bubbleFindingIterations = toCopy->bubbleFindingIterations;
    params->bubbleMinBinomialStrandLikelihood = toCopy->bubbleMinBinomialStrandLikelihood;
    params->bubbleMinBinomialReadSplitLikelihood = toCopy->bubbleMinBinomialReadSplitLikelihood;

    return params;
}

void stRPHmmParameters_parseParametersFromJson(stRPHmmParameters *params, char *buf, size_t r) {
    // Setup parser
    jsmntok_t *tokens;
    char *js;
    int64_t tokenNumber = stJson_setupParser(buf, r, &tokens, &js);

    // Parse tokens, starting at token 1
    // (token 0 is entire object)
    for (int64_t i = 1; i < tokenNumber; i++) {
        jsmntok_t key = tokens[i];
        char *keyString = stJson_token_tostr(js, &key);

        if (strcmp(keyString, "maxNotSumTransitions") == 0) {
            params->maxNotSumTransitions = stJson_parseBool(js, tokens, ++i);
        } else if (strcmp(keyString, "minPartitionsInAColumn") == 0) {
            params->minPartitionsInAColumn = stJson_parseInt(js, tokens, ++i);
        } else if (strcmp(keyString, "maxPartitionsInAColumn") == 0) {
            params->maxPartitionsInAColumn = stJson_parseInt(js, tokens, ++i);
        } else if (strcmp(keyString, "minPosteriorProbabilityForPartition") == 0) {
            params->minPosteriorProbabilityForPartition = stJson_parseFloat(js, tokens, ++i);
        } else if (strcmp(keyString, "maxCoverageDepth") == 0) {
            params->maxCoverageDepth = stJson_parseInt(js, tokens, ++i);
        } else if (strcmp(keyString, "minReadCoverageToSupportPhasingBetweenHeterozygousSites") == 0) {
            params->minReadCoverageToSupportPhasingBetweenHeterozygousSites = stJson_parseInt(js, tokens, ++i);
        } else if (strcmp(keyString, "includeInvertedPartitions") == 0) {
            params->includeInvertedPartitions = stJson_parseBool(js, tokens, ++i);
        } else if (strcmp(keyString, "roundsOfIterativeRefinement") == 0) {
            params->roundsOfIterativeRefinement = stJson_parseInt(js, tokens, ++i);
        } else if (strcmp(keyString, "minPhredScoreForHaplotypePartition") == 0) {
            params->minPhredScoreForHaplotypePartition = stJson_parseInt(js, tokens, ++i);
        } else if (strcmp(keyString, "stitchWithPrimaryReadsOnly") == 0) {
            params->stitchWithPrimaryReadsOnly = stJson_parseBool(js, tokens, ++i);
        } else if (strcmp(keyString, "includeHomozygousVCFEntries") == 0) {
            params->includeHomozygousVCFEntries = stJson_parseBool(js, tokens, ++i);
        } else if (strcmp(keyString, "onlyUsePassVCFEntries") == 0) {
            params->onlyUsePassVCFEntries = stJson_parseBool(js, tokens, ++i);
        } else if (strcmp(keyString, "onlyUseSNPVCFEntries") == 0) {
            params->onlyUseSNPVCFEntries = stJson_parseBool(js, tokens, ++i);
        } else if (strcmp(keyString, "useVariantSelectionAdaptiveSampling") == 0) {
            params->useVariantSelectionAdaptiveSampling = stJson_parseBool(js, tokens, ++i);
        } else if (strcmp(keyString, "variantSelectionAdaptiveSamplingPrimaryThreshold") == 0) {
            params->variantSelectionAdaptiveSamplingPrimaryThreshold = stJson_parseFloat(js, tokens, ++i);
        } else if (strcmp(keyString, "variantSelectionAdaptiveSamplingDesiredBasepairsPerVariant") == 0) {
            params->variantSelectionAdaptiveSamplingDesiredBasepairsPerVariant = stJson_parseInt(js, tokens, ++i);
        } else if (strcmp(keyString, "minVariantQuality") == 0) {
            params->minVariantQuality = stJson_parseFloat(js, tokens, ++i);
        } else if (strcmp(keyString, "updateAllOutputVCFFormatFields") == 0) {
            params->updateAllOutputVCFFormatFields = stJson_parseBool(js, tokens, ++i);
        } else if (strcmp(keyString, "phasesetMinBinomialReadSplitLikelihood") == 0) {
            params->phasesetMinBinomialReadSplitLikelihood = stJson_parseFloat(js, tokens, ++i);
        } else if (strcmp(keyString, "phasesetMaxDiscordantRatio") == 0) {
            params->phasesetMaxDiscordantRatio = stJson_parseFloat(js, tokens, ++i);
        } else if (strcmp(keyString, "bubbleFindingIterations") == 0) {
            params->bubbleFindingIterations = stJson_parseInt(js, tokens, ++i);
        } else if (strcmp(keyString, "bubbleMinBinomialStrandLikelihood") == 0) {
            params->bubbleMinBinomialStrandLikelihood = stJson_parseFloat(js, tokens, ++i);
        } else if (strcmp(keyString, "bubbleMinBinomialReadSplitLikelihood") == 0) {
            params->bubbleMinBinomialReadSplitLikelihood = stJson_parseFloat(js, tokens, ++i);
        } else {
            st_errAbort("ERROR: Unrecognised key in params file: %s\n", keyString);
        }
    }

    // Cleanup
    free(js);
    free(tokens);
}

void stRPHmmParameters_finishParsing(stRPHmmParameters *params) {
    // This is a stub, which is there to check parameters after parsing them out of file(s)
}

/*
 * Params object for polisher
 */

int64_t repeatSubMatrix_parseLogProbabilities(RepeatSubMatrix *repeatSubMatrix, Symbol base, bool strand, char *js,
                                              jsmntok_t *tokens, int64_t tokenIndex) {
    int64_t maxRepeatCount = repeatSubMatrix->maximumRepeatLength;
    int64_t i = stJson_parseFloatArray(repeatSubMatrix_setLogProb(repeatSubMatrix, base, strand, 0, 0),
                                       maxRepeatCount * maxRepeatCount, js, tokens, tokenIndex);
    return i;
}

void repeatSubMatrix_jsonParse(RepeatSubMatrix *repeatSubMatrix, char *buf, size_t r) {
    // Setup parser
    jsmntok_t *tokens;
    char *js;
    int64_t tokenNumber = stJson_setupParser(buf, r, &tokens, &js);

    // TODO: ADD support for Ns
    for (int64_t tokenIndex = 1; tokenIndex < tokenNumber; tokenIndex++) {
        jsmntok_t key = tokens[tokenIndex];
        char *keyString = stJson_token_tostr(js, &key);


        if (stString_eq("baseLogRepeatCounts_AT", keyString)) {
            tokenIndex = stJson_parseFloatArray(repeatSubMatrix->baseLogProbs_AT,
                                                repeatSubMatrix->maximumRepeatLength, js, tokens, tokenIndex + 1);
            continue;
        }

        if (stString_eq("baseLogRepeatCounts_GC", keyString)) {
            tokenIndex = stJson_parseFloatArray(repeatSubMatrix->baseLogProbs_GC,
                                                repeatSubMatrix->maximumRepeatLength, js, tokens, tokenIndex + 1);
            continue;
        }

        if (strlen(keyString) != 31) {
            st_errAbort("ERROR: Unrecognised key in repeat sub matrix json: %s\n", keyString);
        }

        char base = keyString[28];
        if (base != 'A' && base != 'C' && base != 'G' && base != 'T') {
            st_errAbort("ERROR: Unrecognised base in repeat sub matrix json: %s, base=%c\n", keyString, base);
        }
        if (keyString[30] != 'F') {
            st_errAbort("ERROR: Unrecognised strand in repeat sub matrix json: %s, strand:%c\n", keyString,
                        keyString[30]);
        }
        // This sets the probs for the forward strand
        tokenIndex = repeatSubMatrix_parseLogProbabilities(repeatSubMatrix,
                                                           repeatSubMatrix->alphabet->convertCharToSymbol(base), 1, js,
                                                           tokens, tokenIndex + 1);
    }

    // Cleanup
    free(js);
    free(tokens);
}

PolishParams  *polishParams_constructEmpty() {
    // Make empty params object
    PolishParams *params = st_calloc(1, sizeof(PolishParams));

    // Intelligent defaults
    params->useRunLengthEncoding = 1;
    params->referenceBasePenalty = 0.5;
    params->minPosteriorProbForAlignmentAnchors = st_calloc(2, sizeof(double));
    params->minPosteriorProbForAlignmentAnchors[0] = 0.9;
    params->minPosteriorProbForAlignmentAnchors[0] = 10;
    params->minPosteriorProbForAlignmentAnchorsLength = 2;
    params->includeSoftClipping = FALSE;
    params->shuffleChunks = TRUE;
    params->shuffleChunksMethod = SCM_SIZE_DESC;
    params->useRepeatCountsInAlignment = FALSE;
    params->chunkSize = 10000;
    params->chunkBoundary = 1000;
    params->maxDepth = 64;
    params->excessiveDepthThreshold = 512;
    params->includeSecondaryAlignments = FALSE;
    params->includeSupplementaryAlignments = FALSE;
    params->filterAlignmentsWithMapQBelowThisThreshold = 10;
    params->candidateVariantWeight = 0.2;
    params->columnAnchorTrim = 5;
    params->maxConsensusStrings = 100;
    params->repeatSubMatrix = NULL;
    params->stateMachineForGenomeComparison = stateMachine3_constructNucleotide(threeStateAsymmetric);
    params->useReadAlleles = 1;
    params->useReadAllelesInPhasing = 0;
    params->hetSubstitutionProbability = 0.0001;
    params->hetRunLengthSubstitutionProbability = 0.0001;
    params->p = pairwiseAlignmentBandingParameters_construct();

    // At this point the repeat matrix, the hmms for read alignment, the alphabet and the pairwise alignment parameter will be null.

    return params;
}

void polishParams_jsonParse(PolishParams *params, char *buf, size_t r) {
    // Setup parser
    jsmntok_t *tokens;
    char *js;
    int64_t tokenNumber = stJson_setupParser(buf, r, &tokens, &js);

    // Parse tokens, starting at token 1
    // (token 0 is entire object)
    for (int64_t tokenIndex = 1; tokenIndex < tokenNumber; tokenIndex++) {
        jsmntok_t key = tokens[tokenIndex];
        char *keyString = stJson_token_tostr(js, &key);

        if (strcmp(keyString, "useRunLengthEncoding") == 0) {
            params->useRunLengthEncoding = stJson_parseBool(js, tokens, ++tokenIndex);
        } else if (strcmp(keyString, "referenceBasePenalty") == 0) {
            params->referenceBasePenalty = stJson_parseFloat(js, tokens, ++tokenIndex);
        } else if (strcmp(keyString, "minPosteriorProbForAlignmentAnchors") == 0) {
            free(params->minPosteriorProbForAlignmentAnchors); //Cleanup the old one
            int64_t arraySize = tokens[tokenIndex + 1].size;
            if (arraySize % 2 != 0 && arraySize > 0) {
                st_errAbort(
                        "ERROR: length of minPosteriorProbForAlignmentAnchors must be even and greater than zero\n");
            }
            params->minPosteriorProbForAlignmentAnchors = st_calloc(arraySize, sizeof(double));
            params->minPosteriorProbForAlignmentAnchorsLength = arraySize;
            tokenIndex = stJson_parseFloatArray(params->minPosteriorProbForAlignmentAnchors,
                                                arraySize, js, tokens, ++tokenIndex);
            double pValue = 0.0;
            for (int64_t i = 0; i < params->minPosteriorProbForAlignmentAnchorsLength; i += 2) {
                if (params->minPosteriorProbForAlignmentAnchors[i] < pValue ||
                    params->minPosteriorProbForAlignmentAnchors[i] > 1) {
                    st_errAbort(
                            "ERROR: minPosteriorProbForAlignmentAnchors must be even, greater than zero and increasing\n");
                }
                pValue = params->minPosteriorProbForAlignmentAnchors[i];
                if (params->minPosteriorProbForAlignmentAnchors[i + 1] < 0 ||
                    ((int64_t) params->minPosteriorProbForAlignmentAnchors[i + 1]) % 2 != 0) {
                    st_errAbort(
                            "ERROR: minPosteriorProbForAlignmentAnchors diagonal expansion must be, greater than zero and even\n");
                }
            }
        } else if (strcmp(keyString, "repeatCountSubstitutionMatrix") == 0) {
            jsmntok_t tok = tokens[tokenIndex + 1];
            char *tokStr = stJson_token_tostr(js, &tok);
            // TODO: Generalize to not be nucleotide only
            if(params->repeatSubMatrix == NULL) {
                params->repeatSubMatrix = repeatSubMatrix_constructEmpty(alphabet_constructNucleotide());
            }
            repeatSubMatrix_jsonParse(params->repeatSubMatrix, tokStr, strlen(tokStr));
            tokenIndex += stJson_getNestedTokenCount(tokens, tokenIndex + 1);
        } else if (strcmp(keyString, "poaConstructCompareRepeatCounts") == 0) {
            params->poaConstructCompareRepeatCounts = stJson_parseBool(js, tokens, ++tokenIndex);
        } else if (strcmp(keyString, "hmmForwardStrandReadGivenReference") == 0) {
            jsmntok_t tok = tokens[tokenIndex + 1];
            char *tokStr = stJson_token_tostr(js, &tok);
            Hmm *hmmForwardStrandReadGivenReference = hmm_jsonParse(tokStr, strlen(tokStr));
            tokenIndex += stJson_getNestedTokenCount(tokens, tokenIndex + 1);
            // Cleanup any older state machine
            if(params->stateMachineForForwardStrandRead != NULL) {
                stateMachine_destruct(params->stateMachineForForwardStrandRead);
                assert(params->stateMachineForReverseStrandRead != NULL);
                stateMachine_destruct(params->stateMachineForReverseStrandRead);
            }
            params->stateMachineForForwardStrandRead = hmm_getStateMachine(hmmForwardStrandReadGivenReference);
            params->stateMachineForReverseStrandRead = hmm_getStateMachine(hmmForwardStrandReadGivenReference);
            nucleotideEmissions_reverseComplement(
                    (NucleotideEmissions *) params->stateMachineForReverseStrandRead->emissions);
            hmm_destruct(hmmForwardStrandReadGivenReference);
        } else if (strcmp(keyString, "pairwiseAlignmentParameters") == 0) {
            jsmntok_t tok = tokens[tokenIndex + 1];
            char *tokStr = stJson_token_tostr(js, &tok);
            pairwiseAlignmentParameters_jsonParse(params->p, tokStr, strlen(tokStr));
            if (params->p->diagonalExpansion % 2 != 0) {
                st_errAbort("ERROR: pairwiseAlignmentParameters.diagonalExpansion must be even\n");
            }
            tokenIndex += stJson_getNestedTokenCount(tokens, tokenIndex + 1);
        } else if (strcmp(keyString, "shuffleChunks") == 0) {
            params->shuffleChunks = stJson_parseBool(js, tokens, ++tokenIndex);
        } else if (strcmp(keyString, "shuffleChunksMethod") == 0) {
            jsmntok_t tok = tokens[tokenIndex + 1];
            char *tokStr = stJson_token_tostr(js, &tok);
            if (stString_eqcase(tokStr, "random")) {
                params->shuffleChunksMethod = SCM_RANDOM;
            } else if (stString_eqcase(tokStr, "size_desc")) {
                params->shuffleChunksMethod = SCM_SIZE_DESC;
            } else {
                st_errAbort("Invald 'shuffleChunksMethod' parameter '%s'.  Expected ('random', 'size_desc').", tokStr);
            }
            tokenIndex += stJson_getNestedTokenCount(tokens, tokenIndex + 1);
        } else if (strcmp(keyString, "includeSoftClipping") == 0) {
            params->includeSoftClipping = stJson_parseBool(js, tokens, ++tokenIndex);
        } else if (strcmp(keyString, "useRepeatCountsInAlignment") == 0) {
            params->useRepeatCountsInAlignment = stJson_parseBool(js, tokens, ++tokenIndex);
        } else if (strcmp(keyString, "chunkSize") == 0) {
            if (stJson_parseInt(js, tokens, ++tokenIndex) < 0) {
                st_errAbort("ERROR: chunkSize parameter must zero or greater\n");
            }
            params->chunkSize = (uint64_t) stJson_parseInt(js, tokens, tokenIndex);
        } else if (strcmp(keyString, "chunkBoundary") == 0) {
            if (stJson_parseInt(js, tokens, ++tokenIndex) < 0) {
                st_errAbort("ERROR: chunkBoundary parameter must zero or greater\n");
            }
            params->chunkBoundary = (uint64_t) stJson_parseInt(js, tokens, tokenIndex);
        } else if (strcmp(keyString, "maxDepth") == 0) {
            if (stJson_parseInt(js, tokens, ++tokenIndex) < 0) {
                st_errAbort("ERROR: maxDepth parameter must zero or greater\n");
            }
            params->maxDepth = (uint64_t) stJson_parseInt(js, tokens, tokenIndex);
        } else if (strcmp(keyString, "excessiveDepthThreshold") == 0) {
            if (stJson_parseInt(js, tokens, ++tokenIndex) < 0) {
                st_errAbort("ERROR: excessiveDepthThreshold parameter must zero or greater\n");
            }
            params->excessiveDepthThreshold = (uint64_t) stJson_parseInt(js, tokens, tokenIndex);
        } else if (strcmp(keyString, "includeSecondaryAlignments") == 0) {
            params->includeSecondaryAlignments = stJson_parseBool(js, tokens, ++tokenIndex);
        } else if (strcmp(keyString, "includeSupplementaryAlignments") == 0) {
            params->includeSupplementaryAlignments = stJson_parseBool(js, tokens, ++tokenIndex);
        } else if (strcmp(keyString, "filterAlignmentsWithMapQBelowThisThreshold") == 0) {
            if (stJson_parseInt(js, tokens, ++tokenIndex) < 0) {
                st_errAbort("ERROR: filterAlignmentsWithMapQBelowThisThreshold parameter must zero or greater\n");
            }
            params->filterAlignmentsWithMapQBelowThisThreshold = (uint64_t) stJson_parseInt(js, tokens, tokenIndex);
        } else if (strcmp(keyString, "candidateVariantWeight") == 0) {
            if (stJson_parseFloat(js, tokens, ++tokenIndex) < 0) {
                st_errAbort("ERROR: candidateVariantWeight parameter must zero or greater\n");
            }
            params->candidateVariantWeight = stJson_parseFloat(js, tokens, tokenIndex);
        } else if (strcmp(keyString, "columnAnchorTrim") == 0) {
            if (stJson_parseInt(js, tokens, ++tokenIndex) < 0) {
                st_errAbort("ERROR: columnAnchorTrim parameter must zero or greater\n");
            }
            params->columnAnchorTrim = (uint64_t) stJson_parseInt(js, tokens, tokenIndex);
        } else if (strcmp(keyString, "maxConsensusStrings") == 0) {
            if (stJson_parseInt(js, tokens, ++tokenIndex) < 0) {
                st_errAbort("ERROR: maxConsensusStrings parameter must zero or greater\n");
            }
            params->maxConsensusStrings = (uint64_t) stJson_parseInt(js, tokens, tokenIndex);
        } else if (strcmp(keyString, "maxPoaConsensusIterations") == 0) {
            if (stJson_parseInt(js, tokens, ++tokenIndex) < 0) {
                st_errAbort("ERROR: maxPoaConsensusIterations parameter must zero or greater\n");
            }
            params->maxPoaConsensusIterations = (uint64_t) stJson_parseInt(js, tokens, tokenIndex);
        } else if (strcmp(keyString, "minPoaConsensusIterations") == 0) {
            if (stJson_parseInt(js, tokens, ++tokenIndex) < 0) {
                st_errAbort("ERROR: minPoaConsensusIterations parameter must zero or greater\n");
            }
            params->minPoaConsensusIterations = (uint64_t) stJson_parseInt(js, tokens, tokenIndex);
        } else if (strcmp(keyString, "maxRealignmentPolishIterations") == 0) {
            if (stJson_parseInt(js, tokens, ++tokenIndex) < 0) {
                st_errAbort("ERROR: maxRealignmentPolishIterations parameter must zero or greater\n");
            }
            params->maxRealignmentPolishIterations = (uint64_t) stJson_parseInt(js, tokens, tokenIndex);
        } else if (strcmp(keyString, "minRealignmentPolishIterations") == 0) {
            if (stJson_parseInt(js, tokens, ++tokenIndex) < 0) {
                st_errAbort("ERROR: minRealignmentPolishIterations parameter must zero or greater\n");
            }
            params->minRealignmentPolishIterations = (uint64_t) stJson_parseInt(js, tokens, tokenIndex);
        } else if (strcmp(keyString, "filterReadsWhileHaveAtLeastThisCoverage") == 0) {
            if (stJson_parseInt(js, tokens, ++tokenIndex) < 0) {
                st_errAbort("ERROR: filterReadsWhileHaveAtLeastThisCoverage parameter must zero or greater\n");
            }
            params->filterReadsWhileHaveAtLeastThisCoverage = (uint64_t) stJson_parseInt(js, tokens, tokenIndex);
        } else if (strcmp(keyString, "minAvgBaseQuality") == 0) {
            if (stJson_parseFloat(js, tokens, ++tokenIndex) < 0) {
                st_errAbort("ERROR: minAvgBaseQuality parameter must zero or greater\n");
            }
            params->minAvgBaseQuality = stJson_parseFloat(js, tokens, tokenIndex);
        } else if (strcmp(keyString, "hetSubstitutionProbability") == 0) {
            if (stJson_parseFloat(js, tokens, ++tokenIndex) < 0) {
                st_errAbort("ERROR: hetScalingParameter parameter must zero or greater\n");
            }
            params->hetSubstitutionProbability = stJson_parseFloat(js, tokens, tokenIndex);
        } else if (strcmp(keyString, "hetRunLengthSubstitutionProbability") == 0) {
            if (stJson_parseFloat(js, tokens, ++tokenIndex) < 0) {
                st_errAbort("ERROR: hetRunLengthSubstitutionProbability parameter must zero or greater\n");
            }
            params->hetRunLengthSubstitutionProbability = stJson_parseFloat(js, tokens, tokenIndex);
        } else if (strcmp(keyString, "useReadAlleles") == 0) {
            params->useReadAlleles = stJson_parseBool(js, tokens, ++tokenIndex);
        } else if (strcmp(keyString, "skipHaploidPolishingIfDiploid") == 0) {
            params->skipHaploidPolishingIfDiploid = stJson_parseBool(js, tokens, ++tokenIndex);
        } else if (strcmp(keyString, "useReadAllelesInPhasing") == 0) {
            params->useReadAllelesInPhasing = stJson_parseBool(js, tokens, ++tokenIndex);
        } else if (strcmp(keyString, "alphabet") == 0) {
            jsmntok_t tok = tokens[++tokenIndex];
            char *tokStr = stJson_token_tostr(js, &tok);
            if (stString_eq(tokStr, "nucleotide")) {
                params->alphabet = alphabet_constructNucleotide();
            } else {
                st_errAbort("ERROR: Unrecognised alphabet type json: %s\n", tokStr);
            }
        } else {
            st_errAbort("ERROR: Unrecognised key in polish params json: %s\n", keyString);
        }
    }

    // Cleanup
    free(js);
    free(tokens);
}

void polishParams_finishParsing(PolishParams *params) {
    if (params->repeatSubMatrix == NULL && params->useRunLengthEncoding) {
        st_logCritical(
                "  ERROR: Did not find repeat counts specified in json polish params! Will default to MODE estimation\n");
    }
    if (params->stateMachineForForwardStrandRead == NULL) {
        st_errAbort("ERROR: Did not find HMM for alignment of read to a reference specified in json polish params\n");
    }
    if (params->p == NULL) {
        st_errAbort("ERROR: Did not find pairwise alignment params specified in json polish params\n");
    }
    if (params->alphabet == NULL) {
        st_errAbort("ERROR: Did not find alphabet params specified in json polish params\n");
    }

    if (params->useRepeatCountsInAlignment) {
        if (!params->useRunLengthEncoding) {
            st_errAbort(
                    "ERROR: Trying to use repeat counts in read to reference alignment but not using run length encoding\n");
        }
        if (params->repeatSubMatrix == NULL) {
            st_errAbort(
                    "ERROR: Did not find find repeat counts specified in json but trying to use repeat counts in read to reference alignment\n");
        }
        // Replace the emission models of the read to reference state machines
        params->stateMachineForForwardStrandRead->emissions = rleNucleotideEmissions_construct(
                params->stateMachineForForwardStrandRead->emissions, params->repeatSubMatrix, 1);
        params->stateMachineForReverseStrandRead->emissions = rleNucleotideEmissions_construct(
                params->stateMachineForReverseStrandRead->emissions, params->repeatSubMatrix, 0);
    }
}

void polishParams_printParameters(PolishParams *polishParams, FILE *fh) {
    //TODO - complete this - currently this is quite incomplete

    fprintf(fh, "State machine for forward strand read:\n");
    polishParams->stateMachineForForwardStrandRead->printFn(polishParams->stateMachineForForwardStrandRead, fh);

    fprintf(fh, "State machine for reverse strand read:\n");
    polishParams->stateMachineForReverseStrandRead->printFn(polishParams->stateMachineForReverseStrandRead, fh);

    fprintf(fh, "State machine for haplotype comparison:\n");
    polishParams->stateMachineForGenomeComparison->printFn(polishParams->stateMachineForGenomeComparison, fh);
}

void polishParams_destruct(PolishParams *params) {
    if (params->repeatSubMatrix != NULL) repeatSubMatrix_destruct(params->repeatSubMatrix);
    stateMachine_destruct(params->stateMachineForGenomeComparison);
    stateMachine_destruct(params->stateMachineForForwardStrandRead);
    stateMachine_destruct(params->stateMachineForReverseStrandRead);
    pairwiseAlignmentBandingParameters_destruct(params->p);
    free(params->minPosteriorProbForAlignmentAnchors);
    alphabet_destruct(params->alphabet);
    free(params);
}

/*
 * Global params objects
 */

Params *params_constructEmpty() {
    // Make empty params object
    Params *params = st_calloc(1, sizeof(Params));
    params->polishParams = polishParams_constructEmpty();
    params->phaseParams = stRPHmmParameters_construct();

    return params;
}

void params_readParams2(Params *params, char *paramsFile);

void params_jsonParse(Params *params, char *buf, size_t r, char *paramsFile) {
    // Setup parser
    jsmntok_t *tokens;
    char *js;
    int64_t tokenNumber = stJson_setupParser(buf, r, &tokens, &js);

    // Parse tokens, starting at token 1
    // (token 0 is entire object)
    for (int64_t tokenIndex = 1; tokenIndex < tokenNumber; tokenIndex++) {
        jsmntok_t key = tokens[tokenIndex];
        char *keyString = stJson_token_tostr(js, &key);
        if(strcmp(keyString, "include") == 0) {
            jsmntok_t tok = tokens[++tokenIndex];
            char *nestedParamsFile= stJson_token_tostr(js, &tok);
            char *nonRelativeNestedParamsFile = NULL;
            if (nestedParamsFile[0] == '/') {
                nonRelativeNestedParamsFile = stString_copy(nestedParamsFile);
            } else {
                // Make the complete path to the nested params file
                stList *paramsPath = stString_splitByString(paramsFile, "/");
                stList *nestedParamsPath = stString_splitByString(nestedParamsFile, "/");
                free(stList_pop(paramsPath));
                stList_appendAll(paramsPath, nestedParamsPath);
                nonRelativeNestedParamsFile = stString_join2("/", paramsPath);
                st_logDebug(
                        "Parsing nested params file %s, joining params file path: %s and nested params file path: %s\n",
                        nonRelativeNestedParamsFile, paramsFile, nestedParamsFile);
                stList_setDestructor(nestedParamsPath, NULL);
                stList_destruct(nestedParamsPath);
                stList_destruct(paramsPath);
            }
            // Parse the nested params file
            params_readParams2(params, nonRelativeNestedParamsFile);
            //cleanup
            free(nonRelativeNestedParamsFile);
        } else if (strcmp(keyString, "polish") == 0) {
            jsmntok_t tok = tokens[tokenIndex + 1];
            char *tokStr = stJson_token_tostr(js, &tok);
            polishParams_jsonParse(params->polishParams, tokStr, strlen(tokStr));
            tokenIndex += stJson_getNestedTokenCount(tokens, tokenIndex + 1);
        } else if (strcmp(keyString, "phase") == 0) {
            jsmntok_t tok = tokens[tokenIndex + 1];
            char *tokStr = stJson_token_tostr(js, &tok);
            stRPHmmParameters_parseParametersFromJson(params->phaseParams, tokStr, strlen(tokStr));
            tokenIndex += stJson_getNestedTokenCount(tokens, tokenIndex + 1);
        } else {
            st_errAbort("ERROR: Unrecognised key in params json: %s\n", keyString);
        }
    }

    // Cleanup
    free(js);
    free(tokens);
}

void params_readParams2(Params *params, char *paramsFile) {
    // open file and check for existence
    FILE *fh = safe_fopen(paramsFile, "rb");
    if (fh == NULL) {
        st_errAbort("ERROR: Cannot open parameters file %s\n", paramsFile);
    }

    // get file length
    struct stat st;
    fstat(fileno(fh), &st);

    // read file and parse json
    char *buf = st_calloc(st.st_size + 1, sizeof(char));
    size_t readLen = fread(buf, sizeof(char), st.st_size, fh);
    buf[st.st_size] = '\0';
    params_jsonParse(params, buf, readLen, paramsFile);

    // close
    free(buf);
    fclose(fh);
}

Params *params_readParams(char *paramsFile) {
    Params *params = params_constructEmpty();
    params_readParams2(params, paramsFile);
    stRPHmmParameters_finishParsing(params->phaseParams);
    polishParams_finishParsing(params->polishParams); // This initializes everything

    return params;
}

void params_destruct(Params *params) {
    if (params->phaseParams != NULL) stRPHmmParameters_destruct(params->phaseParams);
    if (params->polishParams != NULL) polishParams_destruct(params->polishParams);
    free(params);
}

void params_printParameters(Params *params, FILE *fh) {
    fprintf(fh, "Polish parameters:\n");
    polishParams_printParameters(params->polishParams, fh);
    fprintf(fh, "Phase parameters:\n");
    stRPHmmParameters_printParameters(params->phaseParams, fh);
}
