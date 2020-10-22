/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * pairwiseAligner.h
 *
 *  Created on: 5 Jul 2010
 *      Author: benedictpaten
 */

#ifndef PAIRWISEALIGNER_H_
#define PAIRWISEALIGNER_H_

#include "bioioC.h"
#include "sonLib.h"
#include "pairwiseAlignment.h"
#include "stateMachine.h"

//The exception string
extern const char *PAIRWISE_ALIGNMENT_EXCEPTION_ID;

//Constant that gives the integer value equal to probability 1. Integer probability zero is always 0.
#define PAIR_ALIGNMENT_PROB_1 10000000

typedef struct _pairwiseAlignmentBandingParameters {
    double threshold; //Minimum posterior probability of a match to be added to the output
    int64_t minDiagsBetweenTraceBack; //Minimum x+y diagonals to leave between doing traceback.
    int64_t traceBackDiagonals; //Number of diagonals to leave between trace back diagonal
    int64_t diagonalExpansion; //The number of x-y diagonals to expand around an anchor point
    int64_t constraintDiagonalTrim; //Amount to remove from a diagonal to be considered for a banding constraint
    int64_t splitMatrixBiggerThanThis; //Any matrix in the anchors bigger than this is split into two.
    bool alignAmbiguityCharacters;
    float gapGamma; //The AMAP gap-gamma parameter which controls the degree to which indel probabilities are factored into the alignment.
    bool dynamicAnchorExpansion; // For each alignment anchor specify the expansion of the band individually, instead of using a
    // single expansion
} PairwiseAlignmentParameters;

PairwiseAlignmentParameters *pairwiseAlignmentBandingParameters_construct();

void pairwiseAlignmentBandingParameters_destruct(PairwiseAlignmentParameters *p);

/*
 * Parse from json description. Use defaults specified in pairwiseAlignmentBandingParameters_construct()
 * if a parameter is not specified.
 */
void *pairwiseAlignmentParameters_jsonParse(PairwiseAlignmentParameters *params, char *buf, size_t r);

/*
 * Computes for the forward log probability of aligning the two sequences
 */
double
computeForwardProbability(SymbolString seqX, SymbolString seqY, stList *anchorPairs, PairwiseAlignmentParameters *p,
                          StateMachine *sM,
                          bool alignmentHasRaggedLeftEnd, bool alignmentHasRaggedRightEnd);

/*
 * Gets the set of posterior match probabilities under a simple HMM model of alignment for two DNA sequences.
 */
stList *getAlignedPairs(StateMachine *sM, SymbolString string1, SymbolString string2, PairwiseAlignmentParameters *p,
                        bool alignmentHasRaggedLeftEnd, bool alignmentHasRaggedRightEnd);

/*
 * As getAlignedPairs, but also gives insert and delete probabilities.
 * Return value is by initializing the matches, inserts and deletes lists with values.
 */
void
getAlignedPairsWithIndels(StateMachine *sM, SymbolString string1, SymbolString string2, PairwiseAlignmentParameters *p,
                          stList **alignedPairs, stList **gapXPairs, stList **gapYPairs,
                          bool alignmentHasRaggedLeftEnd, bool alignmentHasRaggedRightEnd);

stList *convertPairwiseForwardStrandAlignmentToAnchorPairs(struct PairwiseAlignment *pA, int64_t trim,
                                                           int64_t diagonalExpansion);

stList *getAlignedPairsUsingAnchors(StateMachine *sM, SymbolString sX, SymbolString sY, stList *anchorPairs,
                                    PairwiseAlignmentParameters *p, bool alignmentHasRaggedLeftEnd,
                                    bool alignmentHasRaggedRightEnd);

void getAlignedPairsWithIndelsUsingAnchors(StateMachine *sM, SymbolString sX, SymbolString sY, stList *anchorPairs,
                                           PairwiseAlignmentParameters *p, stList **alignedPairs, stList **gapXPairs,
                                           stList **gapYPairs,
                                           bool alignmentHasRaggedLeftEnd, bool alignmentHasRaggedRightEnd);

/*
 * As filterPairwiseAlignmentToMakePairsOrdered, but does not use the multiple alignment code. Returns
 * a subset of alignedPairs that form a maximal expected accuracy (MEA) alignment, as described in Schwartz and Pachter.
 * Alignment score is set to the final score of the alignment.
 */
stList *getMaximalExpectedAccuracyPairwiseAlignment(stList *alignedPairs, stList *gapXPairs, stList *gapYPairs,
                                                    int64_t seqXLength, int64_t seqYLength, double *alignmentScore,
                                                    PairwiseAlignmentParameters *p);

/*
 * Shifts pairs in an alignment so that inserts are maximally left shifted.
 */
stList *leftShiftAlignment(stList *alignedPairs, SymbolString seqX, SymbolString seqY);

/*
 * Convenience function that aligns two sequences return a left-shift MEA alignment
 */
stList *
getShiftedMEAAlignment(SymbolString seqX, SymbolString seqY, stList *anchorAlignment, PairwiseAlignmentParameters *p,
                       StateMachine *sM,
                       bool alignmentHasRaggedLeftEnd, bool alignmentHasRaggedRightEnd, double *alignmentScore);

/*
 * Expectation calculation functions for EM algorithms.
 */

void getExpectationsUsingAnchors(StateMachine *sM, Hmm *hmmExpectations, SymbolString sX, SymbolString sY,
                                 stList *anchorPairs,
                                 PairwiseAlignmentParameters *p, bool alignmentHasRaggedLeftEnd,
                                 bool alignmentHasRaggedRightEnd);

void getExpectations(StateMachine *sM, Hmm *hmmExpectations, SymbolString sX, SymbolString sY,
                     PairwiseAlignmentParameters *p, bool alignmentHasRaggedLeftEnd, bool alignmentHasRaggedRightEnd);

/*
 * Methods tested and possibly useful elsewhere
 */

////Diagonal

typedef struct _diagonal {
    int64_t xay; //x + y coordinate
    int64_t xmyL; //smallest x - y coordinate
    int64_t xmyR; //largest x - y coordinate
} Diagonal;

Diagonal diagonal_construct(int64_t xay, int64_t xmyL, int64_t xmyR);

int64_t diagonal_getXay(Diagonal diagonal);

int64_t diagonal_getMinXmy(Diagonal diagonal);

int64_t diagonal_getMaxXmy(Diagonal diagonal);

int64_t diagonal_getWidth(Diagonal diagonal);

int64_t diagonal_getXCoordinate(int64_t xay, int64_t xmy);

int64_t diagonal_getYCoordinate(int64_t xay, int64_t xmy);

int64_t diagonal_equals(Diagonal diagonal1, Diagonal diagonal2);

char *diagonal_getString(Diagonal diagonal);

//Band

typedef struct _band Band;

Band *band_construct(stList *anchorPairs, int64_t lX, int64_t lY,
                     int64_t expansion);

void band_destruct(Band *band);

////Band iterator.

typedef struct _bandIterator BandIterator;

BandIterator *bandIterator_construct(Band *band);

void bandIterator_destruct(BandIterator *bandIterator);

BandIterator *bandIterator_clone(BandIterator *bandIterator);

Diagonal bandIterator_getNext(BandIterator *bandIterator);

Diagonal bandIterator_getPrevious(BandIterator *bandIterator);

//Log add

#define LOG_ZERO -INFINITY

double logAdd(double x, double y);

//Cell calculations

void cell_calculateForward(StateMachine *sM, double *current, double *lower, double *middle, double *upper, Symbol cX,
                           Symbol cY, void *extraArgs);

void cell_calculateBackward(StateMachine *sM, double *current, double *lower, double *middle, double *upper, Symbol cX,
                            Symbol cY, void *extraArgs);

double cell_dotProduct(double *cell1, double *cell2, int64_t stateNumber);

double cell_dotProduct2(double *cell1, StateMachine *sM, double (*getStateValue)(StateMachine *, int64_t));

//DpDiagonal

typedef struct _dpDiagonal DpDiagonal;

DpDiagonal *dpDiagonal_construct(Diagonal diagonal, int64_t stateNumber);

DpDiagonal *dpDiagonal_clone(DpDiagonal *diagonal);

bool dpDiagonal_equals(DpDiagonal *diagonal1, DpDiagonal *diagonal2);

void dpDiagonal_destruct(DpDiagonal *dpDiagonal);

double *dpDiagonal_getCell(DpDiagonal *dpDiagonal, int64_t xmy);

double dpDiagonal_dotProduct(DpDiagonal *diagonal1, DpDiagonal *diagonal2);

void dpDiagonal_zeroValues(DpDiagonal *diagonal);

void
dpDiagonal_initialiseValues(DpDiagonal *diagonal, StateMachine *sM, double (*getStateValue)(StateMachine *, int64_t));

//DpMatrix

typedef struct _dpMatrix DpMatrix;

DpMatrix *dpMatrix_construct(int64_t diagonalNumber, int64_t stateNumber);

void dpMatrix_destruct(DpMatrix *dpMatrix);

DpDiagonal *dpMatrix_getDiagonal(DpMatrix *dpMatrix, int64_t xay);

int64_t dpMatrix_getActiveDiagonalNumber(DpMatrix *dpMatrix);

DpDiagonal *dpMatrix_createDiagonal(DpMatrix *dpMatrix, Diagonal diagonal);

void dpMatrix_deleteDiagonal(DpMatrix *dpMatrix, int64_t xay);

//Diagonal calculations

void diagonalCalculationForward(StateMachine *sM, int64_t xay, DpMatrix *dpMatrix, SymbolString sX, SymbolString sY);

void diagonalCalculationBackward(StateMachine *sM, int64_t xay, DpMatrix *dpMatrix, SymbolString sX, SymbolString sY);

double diagonalCalculationTotalProbability(StateMachine *sM, int64_t xay, DpMatrix *forwardDpMatrix,
                                           DpMatrix *backwardDpMatrix,
                                           SymbolString sX, SymbolString sY);

void diagonalCalculationPosteriorMatchProbs(StateMachine *sM, int64_t xay, DpMatrix *forwardDpMatrix,
                                            DpMatrix *backwardDpMatrix,
                                            SymbolString sX, SymbolString sY,
                                            double totalProbability, PairwiseAlignmentParameters *p, void *extraArgs);

//Banded matrix calculation of posterior probs

void getPosteriorProbsWithBanding(StateMachine *sM, stList *anchorPairs, SymbolString sX, SymbolString sY,
                                  PairwiseAlignmentParameters *p, bool alignmentHasRaggedLeftEnd,
                                  bool alignmentHasRaggedRightEnd,
                                  void (*diagonalPosteriorProbFn)(StateMachine *, int64_t, DpMatrix *, DpMatrix *,
                                                                  SymbolString, SymbolString,
                                                                  double, PairwiseAlignmentParameters *, void *),
                                  void *extraArgs);

//Split over large gaps

stList *getSplitPoints(stList *anchorPairs, int64_t lX, int64_t lY,
                       int64_t maxMatrixSize, bool alignmentHasRaggedLeftEnd, bool alignmentHasRaggedRightEnd);

void getPosteriorProbsWithBandingSplittingAlignmentsByLargeGaps(StateMachine *sM, stList *anchorPairs, SymbolString sX,
                                                                SymbolString sY,
                                                                PairwiseAlignmentParameters *p,
                                                                bool alignmentHasRaggedLeftEnd,
                                                                bool alignmentHasRaggedRightEnd,
                                                                void (*diagonalPosteriorProbFn)(StateMachine *, int64_t,
                                                                                                DpMatrix *, DpMatrix *,
                                                                                                SymbolString,
                                                                                                SymbolString,
                                                                                                double,
                                                                                                PairwiseAlignmentParameters *,
                                                                                                void *),
                                                                void (*coordinateCorrectionFn)(), void *extraArgs);

//Calculate posterior probabilities of being aligned to gaps

int64_t *getIndelProbabilities(stList *alignedPairs, int64_t seqLength, bool xIfTrueElseY);

//Does reweighting, destroys input aligned pairs in the process.
stList *reweightAlignedPairs(stList *alignedPairs,
                             int64_t *indelProbsX, int64_t *indelProbsY, double gapGamma);

stList *reweightAlignedPairs2(stList *alignedPairs, int64_t seqLengthX, int64_t seqLengthY, double gapGamma);

/*
 * Functions to score an alignment by identity / or some proxy to it.
 */

/*
 * Gives the average identity of matches in the alignment, treating indels as mismatches.
 */
int64_t getNumberOfMatchingAlignedPairs(SymbolString subSeqX, SymbolString subSeqY, stList *alignedPairs);

/*
 * Gives the average identity of matches in the alignment, treating indels as mismatches.
 */
double scoreByIdentity(SymbolString subSeqX, SymbolString subSeqY, stList *alignedPairs);

/*
 * Gives the average identity of matches in the alignment, ignoring indels.
 */
double scoreByIdentityIgnoringGaps(SymbolString subSeqX, SymbolString subSeqY, stList *alignedPairs);

/*
 * Gives the average posterior match probability per base of the two sequences, treating bases in indels as having 0 match probability.
 */
double scoreByPosteriorProbability(int64_t lX, int64_t lY, stList *alignedPairs);

/*
 * Gives the average posterior match probability per base of the two sequences, ignoring indels.
 */
double scoreByPosteriorProbabilityIgnoringGaps(stList *alignedPairs);

/*
 * Filter to give co-linear alignment.
 */
stList *filterPairwiseAlignmentToMakePairsOrdered(stList *alignedPairs, SymbolString seqX, SymbolString seqY,
                                                  PairwiseAlignmentParameters *p);

/*
 * Get quick and dirty alignment anchors by finding long chain of shared k-mers, mid-point of each k-mer is then an anchor
 * with given diagonal expansion around it.
 */
stList *getKmerAlignmentAnchors(SymbolString seqX, SymbolString seqY, uint64_t anchorExpansion);

/*
 * Parameters for testing
 */
void setPairwiseAlignerKmerSize(int64_t kmerSize);
void setMinOverlapAnchorPairs(int64_t minOverlapAnchorPairs);

#endif /* PAIRWISEALIGNER_H_ */
