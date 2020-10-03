/*
 * pairwiseAligner.c
 *
 *  Created on: 1 Mar 2012
 *      Author: benedictpaten
 */
//This is being included to make popen work!
#define _XOPEN_SOURCE 500

#include "margin.h"

///////////////////////////////////
///////////////////////////////////
//Diagonal
//
//Structure for working with x-y diagonal of dp matrix
///////////////////////////////////
///////////////////////////////////

const char *PAIRWISE_ALIGNMENT_EXCEPTION_ID = "PAIRWISE_ALIGNMENT_EXCEPTION";

Diagonal diagonal_construct(int64_t xay, int64_t xmyL, int64_t xmyR) {
    if ((xay + xmyL) % 2 != 0 || (xay + xmyR) % 2 != 0 || xmyL > xmyR) {
        stThrowNew(PAIRWISE_ALIGNMENT_EXCEPTION_ID,
                   "Attempt to create diagonal with invalid coordinates: xay %" PRIi64 " xmyL %" PRIi64 " xmyR %" PRIi64 "",
                   xay, xmyL, xmyR);
    }
    Diagonal diagonal;
    diagonal.xay = xay;
    diagonal.xmyL = xmyL;
    diagonal.xmyR = xmyR;
    assert(xmyL <= xmyR);
    assert(xay >= 0);
    return diagonal;
}

inline int64_t diagonal_getXay(Diagonal diagonal) {
    return diagonal.xay;
}

inline int64_t diagonal_getMinXmy(Diagonal diagonal) {
    return diagonal.xmyL;
}

inline int64_t diagonal_getMaxXmy(Diagonal diagonal) {
    return diagonal.xmyR;
}

inline int64_t diagonal_getWidth(Diagonal diagonal) {
    return (diagonal.xmyR - diagonal.xmyL) / 2 + 1;
}

inline int64_t diagonal_getXCoordinate(int64_t xay, int64_t xmy) {
    assert((xay + xmy) % 2 == 0);
    return (xay + xmy) / 2;
}

inline int64_t diagonal_equals(Diagonal diagonal1, Diagonal diagonal2) {
    return diagonal1.xay == diagonal2.xay && diagonal1.xmyL == diagonal2.xmyL && diagonal1.xmyR == diagonal2.xmyR;
}

inline int64_t diagonal_getYCoordinate(int64_t xay, int64_t xmy) {
    assert((xay - xmy) % 2 == 0);
    return (xay - xmy) / 2;
}

inline char *diagonal_getString(Diagonal diagonal) {
    return stString_print("Diagonal, xay: %" PRIi64 " xmyL %" PRIi64 ", xmyR: %" PRIi64 "", diagonal_getXay(diagonal),
                          diagonal_getMinXmy(diagonal), diagonal_getMaxXmy(diagonal));
}

///////////////////////////////////
///////////////////////////////////
//Band Iterator
//
//Iterator for walking along x+y diagonals in banded fashion
//(using a set of anchor constraints)
///////////////////////////////////
///////////////////////////////////

struct _band {
    Diagonal *diagonals;
    int64_t lXalY;
};

static int64_t band_avoidOffByOne(int64_t xay, int64_t xmy) {
    return (xay + xmy) % 2 == 0 ? xmy : xmy + 1;
}

static void band_setCurrentDiagonalP(int64_t *xmy, int64_t i, int64_t j, int64_t k) {
    if (i < j) {
        *xmy += (int64_t) (2 * ((int64_t) j - i) * (int64_t) k);
    }
}

static Diagonal band_setCurrentDiagonal(int64_t xay, int64_t xL, int64_t yL, int64_t xU, int64_t yU) {
    int64_t xmyL = xL - yL;
    int64_t xmyR = xU - yU;

    assert(xay >= xL + yU);
    assert(xay <= xU + yL);

    //Avoid in-between undefined x,y coordinate positions when intersecting xay and xmy.
    xmyL = band_avoidOffByOne(xay, xmyL);
    xmyR = band_avoidOffByOne(xay, xmyR);

    //Bound the xmy coordinates by the xL, yL and xU, yU band boundaries
    band_setCurrentDiagonalP(&xmyL, diagonal_getXCoordinate(xay, xmyL), xL, 1);
    band_setCurrentDiagonalP(&xmyL, yL, diagonal_getYCoordinate(xay, xmyL), 1);
    band_setCurrentDiagonalP(&xmyR, xU, diagonal_getXCoordinate(xay, xmyR), -1);
    band_setCurrentDiagonalP(&xmyR, diagonal_getYCoordinate(xay, xmyR), yU, -1);

    return diagonal_construct(xay, xmyL, xmyR);
}

static int64_t band_boundCoordinate(int64_t z, int64_t lZ) {
    return z < 0 ? 0 : (z > lZ ? lZ : z);
}

Band *band_constructDynamic(stList *anchorPairs, int64_t lX, int64_t lY) {
    //Prerequisities
    assert(lX >= 0);
    assert(lY >= 0);

    Band *band = st_malloc(sizeof(Band));
    band->diagonals = st_malloc(sizeof(Diagonal) * (lX + lY + 1));
    band->lXalY = lX + lY;

    //Now initialise the diagonals
    int64_t anchorPairIndex = 0;
    int64_t xay = 0;
    int64_t pxay = 0, pxmy = 0;
    int64_t nxay = 0, nxmy = 0;
    int64_t xL = 0, yL = 0, xU = 0, yU = 0, expansion = 0;

    while (xay <= band->lXalY) {
        band->diagonals[xay] = band_setCurrentDiagonal(xay, xL, yL, xU, yU);
        if (nxay == xay++) {
            //The previous diagonals become the next
            pxay = nxay;
            pxmy = nxmy;

            int64_t x = lX, y = lY;
            if (anchorPairIndex < stList_length(anchorPairs)) {
                stIntTuple *anchorPair = stList_get(anchorPairs, anchorPairIndex++);
                x = stIntTuple_get(anchorPair, 0) + 1; //Plus ones, because matrix coordinates are +1 the sequence ones
                y = stIntTuple_get(anchorPair, 1) + 1;
                expansion = stIntTuple_get(anchorPair, 2);

                //Check the anchor pairs
                assert(x > diagonal_getXCoordinate(pxay, pxmy));
                assert(y > diagonal_getYCoordinate(pxay, pxmy));
                assert(x <= lX);
                assert(y <= lY);
                assert(x > 0);
                assert(y > 0);
                assert(expansion >= 0);
                assert(expansion % 2 == 0);
            }

            nxay = x + y;
            nxmy = x - y;

            //Now call to set the lower and upper x,y coordinates
            xL = band_boundCoordinate(diagonal_getXCoordinate(pxay, pxmy - expansion), lX);
            yL = band_boundCoordinate(diagonal_getYCoordinate(nxay, nxmy - expansion), lY);
            xU = band_boundCoordinate(diagonal_getXCoordinate(nxay, nxmy + expansion), lX);
            yU = band_boundCoordinate(diagonal_getYCoordinate(pxay, pxmy + expansion), lY);
        }
    }

    return band;
}

Band *band_construct(stList *anchorPairs, int64_t lX, int64_t lY, int64_t expansion) {
    //Prerequisities
    assert(lX >= 0);
    assert(lY >= 0);
    assert(expansion % 2 == 0);

    Band *band = st_malloc(sizeof(Band));
    band->diagonals = st_malloc(sizeof(Diagonal) * (lX + lY + 1));
    band->lXalY = lX + lY;

    //Now initialise the diagonals
    int64_t anchorPairIndex = 0;
    int64_t xay = 0;
    int64_t pxay = 0, pxmy = 0;
    int64_t nxay = 0, nxmy = 0;
    int64_t xL = 0, yL = 0, xU = 0, yU = 0;

    while (xay <= band->lXalY) {
        band->diagonals[xay] = band_setCurrentDiagonal(xay, xL, yL, xU, yU);
        if (nxay == xay++) {
            //The previous diagonals become the next
            pxay = nxay;
            pxmy = nxmy;

            int64_t x = lX, y = lY;
            if (anchorPairIndex < stList_length(anchorPairs)) {
                stIntTuple *anchorPair = stList_get(anchorPairs, anchorPairIndex++);
                x = stIntTuple_get(anchorPair, 0) + 1; //Plus ones, because matrix coordinates are +1 the sequence ones
                y = stIntTuple_get(anchorPair, 1) + 1;

                //Check the anchor pairs
                assert(x > diagonal_getXCoordinate(pxay, pxmy));
                assert(y > diagonal_getYCoordinate(pxay, pxmy));
                assert(x <= lX);
                assert(y <= lY);
                assert(x > 0);
                assert(y > 0);
            }

            nxay = x + y;
            nxmy = x - y;

            //Now call to set the lower and upper x,y coordinates
            xL = band_boundCoordinate(diagonal_getXCoordinate(pxay, pxmy - expansion), lX);
            yL = band_boundCoordinate(diagonal_getYCoordinate(nxay, nxmy - expansion), lY);
            xU = band_boundCoordinate(diagonal_getXCoordinate(nxay, nxmy + expansion), lX);
            yU = band_boundCoordinate(diagonal_getYCoordinate(pxay, pxmy + expansion), lY);
        }
    }

    return band;
}

void band_destruct(Band *band) {
    free(band->diagonals);
    free(band);
}

struct _bandIterator {
    Band *band;
    int64_t index;
};

BandIterator *bandIterator_construct(Band *band) {
    BandIterator *bandIterator = st_malloc(sizeof(BandIterator));
    bandIterator->band = band;
    bandIterator->index = 0;
    return bandIterator;
}

BandIterator *bandIterator_clone(BandIterator *bandIterator) {
    BandIterator *bandIterator2 = st_malloc(sizeof(BandIterator));
    memcpy(bandIterator2, bandIterator, sizeof(BandIterator));
    return bandIterator2;
}

void bandIterator_destruct(BandIterator *bandIterator) {
    free(bandIterator);
}

Diagonal bandIterator_getNext(BandIterator *bandIterator) {
    Diagonal diagonal = bandIterator->band->diagonals[
            bandIterator->index > bandIterator->band->lXalY ? bandIterator->band->lXalY : bandIterator->index];
    if (bandIterator->index <= bandIterator->band->lXalY) {
        bandIterator->index++;
    }
    return diagonal;
}

Diagonal bandIterator_getPrevious(BandIterator *bandIterator) {
    if (bandIterator->index > 0) {
        bandIterator->index--;
    }
    return bandIterator->band->diagonals[bandIterator->index];
}

///////////////////////////////////
///////////////////////////////////
//Log Add functions
//
//Interpolation function for doing log add
///////////////////////////////////
///////////////////////////////////

#define logUnderflowThreshold 7.5
#define posteriorMatchThreshold 0.01

static inline double lookup(double x) {
    //return log (exp (x) + 1);
    assert(x >= 0.00f);
    assert(x <= logUnderflowThreshold);
    if (x <= 1.00f)
        return ((-0.009350833524763f * x + 0.130659527668286f) * x + 0.498799810682272f) * x + 0.693203116424741f;
    if (x <= 2.50f)
        return ((-0.014532321752540f * x + 0.139942324101744f) * x + 0.495635523139337f) * x + 0.692140569840976f;
    if (x <= 4.50f)
        return ((-0.004605031767994f * x + 0.063427417320019f) * x + 0.695956496475118f) * x + 0.514272634594009f;
    return ((-0.000458661602210f * x + 0.009695946122598f) * x + 0.930734667215156f) * x + 0.168037164329057f;
}

double logAdd(double x, double y) {
    if (x < y)
        return (x == LOG_ZERO || y - x >= logUnderflowThreshold) ? y : lookup(y - x) + x;
    return (y == LOG_ZERO || x - y >= logUnderflowThreshold) ? x : lookup(x - y) + y;
}

///////////////////////////////////
///////////////////////////////////
//Cell calculations
//
//A cell is a set of states associated with an x, y coordinate.
//These functions do the forward/backward calculations for the pairwise
//alignment model.
///////////////////////////////////
///////////////////////////////////

static inline void doTransitionForward(double *fromCells, double *toCells, int64_t from, int64_t to, double eP,
                                       double tP, void *extraArgs) {
    toCells[to] = logAdd(toCells[to], fromCells[from] + (eP + tP));
}

void cell_calculateForward(StateMachine *sM, double *current, double *lower, double *middle, double *upper, Symbol cX,
                           Symbol cY,
                           void *extraArgs) {
    sM->cellCalculate(sM, current, lower, middle, upper, cX, cY, doTransitionForward, extraArgs);
}

static inline void doTransitionBackward(double *fromCells, double *toCells, int64_t from, int64_t to, double eP,
                                        double tP, void *extraArgs) {
    fromCells[from] = logAdd(fromCells[from], toCells[to] + (eP + tP));
}

void cell_calculateBackward(StateMachine *sM, double *current, double *lower, double *middle, double *upper, Symbol cX,
                            Symbol cY,
                            void *extraArgs) {
    sM->cellCalculate(sM, current, lower, middle, upper, cX, cY, doTransitionBackward, extraArgs);
}

double cell_dotProduct(double *cell1, double *cell2, int64_t stateNumber) {
    double totalProb = cell1[0] + cell2[0];
    for (int64_t i = 1; i < stateNumber; i++) {
        totalProb = logAdd(totalProb, cell1[i] + cell2[i]);
    }
    return totalProb;
}

double cell_dotProduct2(double *cell, StateMachine *sM, double (*getStateValue)(StateMachine *, int64_t)) {
    double totalProb = cell[0] + getStateValue(sM, 0);
    for (int64_t i = 1; i < sM->stateNumber; i++) {
        totalProb = logAdd(totalProb, cell[i] + getStateValue(sM, i));
    }
    return totalProb;
}

static inline void updateExpectations(double *fromCells, double *toCells, int64_t from, int64_t to, double eP,
                                      double tP, void *extraArgs) {
    //void *extraArgs2[2] = { &totalProbability, hmmExpectations };
    double totalProbability = *((double *) ((void **) extraArgs)[0]);
    Hmm *hmmExpectations = ((void **) extraArgs)[1];
    Symbol x = *((Symbol *) ((void **) extraArgs)[2]);
    Symbol y = *((Symbol *) ((void **) extraArgs)[3]);
    //Calculate posterior probability of the transition/emission pair
    double p = exp(fromCells[from] + toCells[to] + (eP + tP) - totalProbability);
    //Add in the expectation of the transition
    hmm_addToTransitionExpectation(hmmExpectations, from, to, p);

    // Training of emissions has been temporarily removed
    //if(x < SYMBOL_NUMBER_NO_N && y < SYMBOL_NUMBER_NO_N) { //Ignore gaps involving Ns.
    //    hmm_addToEmissionsExpectation(hmmExpectations, to, x, y, p);
    //}
}

static void
cell_calculateExpectation(StateMachine *sM, double *current, double *lower, double *middle, double *upper, Symbol cX,
                          Symbol cY,
                          void *extraArgs) {
    void *extraArgs2[4] = {((void **) extraArgs)[0], ((void **) extraArgs)[1], &cX, &cY};
    sM->cellCalculate(sM, current, lower, middle, upper, cX, cY, updateExpectations, extraArgs2);
}

///////////////////////////////////
///////////////////////////////////
//DpDiagonal
//
//Structure for storing a x-y diagonal of the dp matrix
///////////////////////////////////
///////////////////////////////////

struct _dpDiagonal {
    Diagonal diagonal;
    int64_t stateNumber;
    double *cells;
};

DpDiagonal *dpDiagonal_construct(Diagonal diagonal, int64_t stateNumber) {
    DpDiagonal *dpDiagonal = st_malloc(sizeof(DpDiagonal));
    dpDiagonal->diagonal = diagonal;
    dpDiagonal->stateNumber = stateNumber;
    assert(diagonal_getWidth(diagonal) >= 0);
    dpDiagonal->cells = st_malloc(sizeof(double) * stateNumber * (int64_t) diagonal_getWidth(diagonal));
    return dpDiagonal;
}

DpDiagonal *dpDiagonal_clone(DpDiagonal *diagonal) {
    DpDiagonal *diagonal2 = dpDiagonal_construct(diagonal->diagonal, diagonal->stateNumber);
    memcpy(diagonal2->cells, diagonal->cells,
           sizeof(double) * diagonal_getWidth(diagonal->diagonal) * diagonal->stateNumber);
    return diagonal2;
}

bool dpDiagonal_equals(DpDiagonal *diagonal1, DpDiagonal *diagonal2) {
    if (!diagonal_equals(diagonal1->diagonal, diagonal2->diagonal)) {
        return 0;
    }
    if (diagonal1->stateNumber != diagonal2->stateNumber) {
        return 0;
    }
    for (int64_t i = 0; i < diagonal_getWidth(diagonal1->diagonal) * diagonal1->stateNumber; i++) {
        if (diagonal1->cells[i] != diagonal2->cells[i]) {
            return 0;
        }
    }
    return 1;
}

void dpDiagonal_destruct(DpDiagonal *dpDiagonal) {
    free(dpDiagonal->cells);
    free(dpDiagonal);
}

double *dpDiagonal_getCell(DpDiagonal *dpDiagonal, int64_t xmy) {
    if (xmy < dpDiagonal->diagonal.xmyL || xmy > dpDiagonal->diagonal.xmyR) {
        return NULL;
    }
    assert((diagonal_getXay(dpDiagonal->diagonal) + xmy) % 2 == 0);
    return &dpDiagonal->cells[((xmy - dpDiagonal->diagonal.xmyL) / 2) * dpDiagonal->stateNumber];
}

void dpDiagonal_zeroValues(DpDiagonal *diagonal) {
    for (int64_t i = 0; i < diagonal_getWidth(diagonal->diagonal) * diagonal->stateNumber; i++) {
        diagonal->cells[i] = LOG_ZERO;
    }
}

void
dpDiagonal_initialiseValues(DpDiagonal *diagonal, StateMachine *sM, double (*getStateValue)(StateMachine *, int64_t)) {
    for (int64_t i = diagonal_getMinXmy(diagonal->diagonal); i <= diagonal_getMaxXmy(diagonal->diagonal); i += 2) {
        double *cell = dpDiagonal_getCell(diagonal, i);
        assert(cell != NULL);
        for (int64_t j = 0; j < diagonal->stateNumber; j++) {
            cell[j] = getStateValue(sM, j);
        }
    }
}

double dpDiagonal_dotProduct(DpDiagonal *diagonal1, DpDiagonal *diagonal2) {
    double totalProbability = LOG_ZERO;
    Diagonal diagonal = diagonal1->diagonal;
    int64_t xmy = diagonal_getMinXmy(diagonal);
    while (xmy <= diagonal_getMaxXmy(diagonal)) {
        totalProbability = logAdd(totalProbability,
                                  cell_dotProduct(dpDiagonal_getCell(diagonal1, xmy),
                                                  dpDiagonal_getCell(diagonal2, xmy), diagonal1->stateNumber));
        xmy += 2;
    }
    return totalProbability;
}

///////////////////////////////////
///////////////////////////////////
//DpMatrix
//
//Structure for storing dp-matrix
///////////////////////////////////
///////////////////////////////////

struct _dpMatrix {
    DpDiagonal **diagonals;
    int64_t diagonalNumber;
    int64_t activeDiagonals;
    int64_t stateNumber;
};

DpMatrix *dpMatrix_construct(int64_t diagonalNumber, int64_t stateNumber) {
    assert(diagonalNumber >= 0);
    DpMatrix *dpMatrix = st_malloc(sizeof(DpMatrix));
    dpMatrix->diagonalNumber = diagonalNumber;
    dpMatrix->diagonals = st_calloc(dpMatrix->diagonalNumber + 1, sizeof(DpDiagonal *));
    dpMatrix->activeDiagonals = 0;
    dpMatrix->stateNumber = stateNumber;
    return dpMatrix;
}

void dpMatrix_destruct(DpMatrix *dpMatrix) {
    assert(dpMatrix->activeDiagonals == 0);
    free(dpMatrix->diagonals);
    free(dpMatrix);
}

DpDiagonal *dpMatrix_getDiagonal(DpMatrix *dpMatrix, int64_t xay) {
    if (xay < 0 || xay > dpMatrix->diagonalNumber) {
        return NULL;
    }
    return dpMatrix->diagonals[xay];
}

int64_t dpMatrix_getActiveDiagonalNumber(DpMatrix *dpMatrix) {
    return dpMatrix->activeDiagonals;
}

DpDiagonal *dpMatrix_createDiagonal(DpMatrix *dpMatrix, Diagonal diagonal) {
    assert(diagonal.xay >= 0);
    assert(diagonal.xay <= dpMatrix->diagonalNumber);
    assert(dpMatrix_getDiagonal(dpMatrix, diagonal.xay) == NULL);
    DpDiagonal *dpDiagonal = dpDiagonal_construct(diagonal, dpMatrix->stateNumber);
    dpMatrix->diagonals[diagonal_getXay(diagonal)] = dpDiagonal;
    dpMatrix->activeDiagonals++;
    return dpDiagonal;
}

void dpMatrix_deleteDiagonal(DpMatrix *dpMatrix, int64_t xay) {
    assert(xay >= 0);
    assert(xay <= dpMatrix->diagonalNumber);
    if (dpMatrix->diagonals[xay] != NULL) {
        dpMatrix->activeDiagonals--;
        assert(dpMatrix->activeDiagonals >= 0);
        dpDiagonal_destruct(dpMatrix->diagonals[xay]);
        dpMatrix->diagonals[xay] = NULL;
    }
}

///////////////////////////////////
///////////////////////////////////
//Diagonal DP Calculations
//
//Functions which do forward/backward/posterior calculations
//between diagonal rows of a dp-matrix
///////////////////////////////////
///////////////////////////////////

static Symbol getXCharacter(const SymbolString sX, int64_t xay, int64_t xmy) {
    int64_t x = diagonal_getXCoordinate(xay, xmy);
    assert(x >= 0 && x <= sX.length);
    return x > 0 ? sX.sequence[x - 1] : 4; //n; TODO: this is a hack, must fix
}

static Symbol getYCharacter(const SymbolString sY, int64_t xay, int64_t xmy) {
    int64_t y = diagonal_getYCoordinate(xay, xmy);
    assert(y >= 0 && y <= sY.length);
    return y > 0 ? sY.sequence[y - 1] : 4; //n; TODO: this is a hack, must fix
}

static void
diagonalCalculation(StateMachine *sM, DpDiagonal *dpDiagonal, DpDiagonal *dpDiagonalM1, DpDiagonal *dpDiagonalM2,
                    const SymbolString sX, const SymbolString sY,
                    void (*cellCalculation)(StateMachine *, double *, double *, double *, double *, Symbol, Symbol,
                                            void *), void *extraArgs) {
    Diagonal diagonal = dpDiagonal->diagonal;
    int64_t xmy = diagonal_getMinXmy(diagonal);
    while (xmy <= diagonal_getMaxXmy(diagonal)) {
        Symbol x = getXCharacter(sX, diagonal_getXay(diagonal), xmy);
        Symbol y = getYCharacter(sY, diagonal_getXay(diagonal), xmy);
        double *current = dpDiagonal_getCell(dpDiagonal, xmy);
        double *lower = dpDiagonalM1 == NULL ? NULL : dpDiagonal_getCell(dpDiagonalM1, xmy - 1);
        double *middle = dpDiagonalM2 == NULL ? NULL : dpDiagonal_getCell(dpDiagonalM2, xmy);
        double *upper = dpDiagonalM1 == NULL ? NULL : dpDiagonal_getCell(dpDiagonalM1, xmy + 1);
        cellCalculation(sM, current, lower, middle, upper, x, y, extraArgs);
        xmy += 2;
    }
}

void diagonalCalculationForward(StateMachine *sM, int64_t xay, DpMatrix *dpMatrix, const SymbolString sX,
                                const SymbolString sY) {
    diagonalCalculation(sM, dpMatrix_getDiagonal(dpMatrix, xay), dpMatrix_getDiagonal(dpMatrix, xay - 1),
                        dpMatrix_getDiagonal(dpMatrix, xay - 2), sX, sY, cell_calculateForward, NULL);
}

void diagonalCalculationBackward(StateMachine *sM, int64_t xay, DpMatrix *dpMatrix, const SymbolString sX,
                                 const SymbolString sY) {
    diagonalCalculation(sM, dpMatrix_getDiagonal(dpMatrix, xay), dpMatrix_getDiagonal(dpMatrix, xay - 1),
                        dpMatrix_getDiagonal(dpMatrix, xay - 2), sX, sY, cell_calculateBackward, NULL);
}

double diagonalCalculationTotalProbability(StateMachine *sM, int64_t xay, DpMatrix *forwardDpMatrix,
                                           DpMatrix *backwardDpMatrix,
                                           const SymbolString sX, const SymbolString sY) {
    //Get the forward and backward diagonals
    DpDiagonal *forwardDiagonal = dpMatrix_getDiagonal(forwardDpMatrix, xay);
    DpDiagonal *backDiagonal = dpMatrix_getDiagonal(backwardDpMatrix, xay);
    double totalProbability = dpDiagonal_dotProduct(forwardDiagonal, backDiagonal);
    //Now calculate the contribution of matches through xay.
    forwardDiagonal = dpMatrix_getDiagonal(forwardDpMatrix, xay - 1);
    backDiagonal = dpMatrix_getDiagonal(backwardDpMatrix, xay + 1);
    if (backDiagonal != NULL && forwardDiagonal != NULL) {
        DpDiagonal *matchDiagonal = dpDiagonal_clone(backDiagonal);
        dpDiagonal_zeroValues(matchDiagonal);
        diagonalCalculation(sM, matchDiagonal, NULL, forwardDiagonal, sX, sY, cell_calculateForward, NULL);
        totalProbability = logAdd(totalProbability, dpDiagonal_dotProduct(matchDiagonal, backDiagonal));
        dpDiagonal_destruct(matchDiagonal);
    }
    return totalProbability;
}

void addPosteriorProb(int64_t x, int64_t y, double posteriorProbability, stList *posteriorProbs,
                      PairwiseAlignmentParameters *p) {
    if (posteriorProbability >= p->threshold) {
        if (posteriorProbability > 1.0) {
            posteriorProbability = 1.0;
        }
        posteriorProbability = floor(posteriorProbability * PAIR_ALIGNMENT_PROB_1);

        stList_append(posteriorProbs, stIntTuple_construct3((int64_t) posteriorProbability, x - 1, y - 1));
    }
}

void diagonalCalculationPosteriorMatchProbs(StateMachine *sM, int64_t xay, DpMatrix *forwardDpMatrix,
                                            DpMatrix *backwardDpMatrix,
                                            const SymbolString sX, const SymbolString sY, double totalProbability,
                                            PairwiseAlignmentParameters *p,
                                            void *extraArgs) {
    assert(p->threshold >= 0.0);
    assert(p->threshold <= 1.0);
    stList *alignedPairs = ((void **) extraArgs)[0];
    DpDiagonal *forwardDiagonal = dpMatrix_getDiagonal(forwardDpMatrix, xay);
    DpDiagonal *backDiagonal = dpMatrix_getDiagonal(backwardDpMatrix, xay);
    Diagonal diagonal = forwardDiagonal->diagonal;
    int64_t xmy = diagonal_getMinXmy(diagonal);
    //Walk over the cells computing the posteriors
    while (xmy <= diagonal_getMaxXmy(diagonal)) {
        int64_t x = diagonal_getXCoordinate(diagonal_getXay(diagonal), xmy);
        int64_t y = diagonal_getYCoordinate(diagonal_getXay(diagonal), xmy);
        if (x > 0 && y > 0) {
            double *cellForward = dpDiagonal_getCell(forwardDiagonal, xmy);
            double *cellBackward = dpDiagonal_getCell(backDiagonal, xmy);
            double posteriorProbability = exp(
                    (cellForward[sM->matchState] + cellBackward[sM->matchState]) - totalProbability);
            addPosteriorProb(x, y, posteriorProbability, alignedPairs, p);
        }
        xmy += 2;
    }
}

void
diagonalCalculationPosteriorProbs(StateMachine *sM, int64_t xay, DpMatrix *forwardDpMatrix, DpMatrix *backwardDpMatrix,
                                  const SymbolString sX, const SymbolString sY, double totalProbability,
                                  PairwiseAlignmentParameters *p,
                                  void *extraArgs) {
    assert(p->threshold >= 0.0);
    assert(p->threshold <= 1.0);

    stList *alignedPairs = ((void **) extraArgs)[0];
    stList *gapXPairs = ((void **) extraArgs)[2];
    stList *gapYPairs = ((void **) extraArgs)[4];

    DpDiagonal *forwardDiagonal = dpMatrix_getDiagonal(forwardDpMatrix, xay);
    DpDiagonal *backDiagonal = dpMatrix_getDiagonal(backwardDpMatrix, xay);
    Diagonal diagonal = forwardDiagonal->diagonal;
    int64_t xmy = diagonal_getMinXmy(diagonal);
    //Walk over the cells computing the posteriors
    while (xmy <= diagonal_getMaxXmy(diagonal)) {
        int64_t x = diagonal_getXCoordinate(diagonal_getXay(diagonal), xmy);
        int64_t y = diagonal_getYCoordinate(diagonal_getXay(diagonal), xmy);

        double *cellForward = dpDiagonal_getCell(forwardDiagonal, xmy);
        double *cellBackward = dpDiagonal_getCell(backDiagonal, xmy);
        if (x > 0 && y > 0) {
            // Posterior match prob
            double posteriorProbability = exp(
                    (cellForward[sM->matchState] + cellBackward[sM->matchState]) - totalProbability);
            addPosteriorProb(x, y, posteriorProbability, alignedPairs, p);
        }

        if (x > 0) {
            double posteriorProbability = exp(
                    (cellForward[sM->gapXState] + cellBackward[sM->gapXState]) - totalProbability);
            addPosteriorProb(x, y, posteriorProbability, gapXPairs, p);
        }

        if (y > 0) {
            double posteriorProbability = exp(
                    (cellForward[sM->gapYState] + cellBackward[sM->gapYState]) - totalProbability);
            addPosteriorProb(x, y, posteriorProbability, gapYPairs, p);
        }

        xmy += 2;
    }
}

static void
diagonalCalculationExpectations(StateMachine *sM, int64_t xay, DpMatrix *forwardDpMatrix, DpMatrix *backwardDpMatrix,
                                const SymbolString sX, const SymbolString sY, double totalProbability,
                                PairwiseAlignmentParameters *p,
                                void *extraArgs) {
    /*
     * Updates the expectations of the transitions/emissions for the given diagonal.
     */
    Hmm *hmmExpectations = extraArgs;
    void *extraArgs2[2] = {&totalProbability, hmmExpectations};
    hmmExpectations->likelihood += totalProbability; //We do this once per diagonal, which is a hack, rather than for the whole matrix. The correction factor is approximately 1/number of diagonals.
    diagonalCalculation(sM, dpMatrix_getDiagonal(backwardDpMatrix, xay), dpMatrix_getDiagonal(forwardDpMatrix, xay - 1),
                        dpMatrix_getDiagonal(forwardDpMatrix, xay - 2), sX, sY, cell_calculateExpectation, extraArgs2);
}

///////////////////////////////////
///////////////////////////////////
//Banded alignment routine to calculate posterior match probs
//
//
///////////////////////////////////
///////////////////////////////////

void getPosteriorProbsWithBanding(StateMachine *sM, stList *anchorPairs, const SymbolString sX, const SymbolString sY,
                                  PairwiseAlignmentParameters *p, bool alignmentHasRaggedLeftEnd,
                                  bool alignmentHasRaggedRightEnd,
                                  void (*diagonalPosteriorProbFn)(StateMachine *, int64_t, DpMatrix *, DpMatrix *,
                                                                  const SymbolString, const SymbolString, double,
                                                                  PairwiseAlignmentParameters *, void *),
                                  void *extraArgs) {
    //Prerequisites
    assert(p->traceBackDiagonals >= 1);
    assert(p->diagonalExpansion >= 0);
    assert(p->diagonalExpansion % 2 == 0);
    assert(p->minDiagsBetweenTraceBack >= 2);
    assert(p->traceBackDiagonals + 1 < p->minDiagsBetweenTraceBack);

    int64_t diagonalNumber = sX.length + sY.length;
    if (diagonalNumber == 0) { //Deal with trivial case
        return;
    }

    //Primitives for the forward matrix recursion
    Band *band = p->dynamicAnchorExpansion ? band_constructDynamic(anchorPairs, sX.length, sY.length) : band_construct(
            anchorPairs, sX.length, sY.length, p->diagonalExpansion);
    BandIterator *forwardBandIterator = bandIterator_construct(band);
    DpMatrix *forwardDpMatrix = dpMatrix_construct(diagonalNumber, sM->stateNumber);
    dpDiagonal_initialiseValues(dpMatrix_createDiagonal(forwardDpMatrix, bandIterator_getNext(forwardBandIterator)), sM,
                                alignmentHasRaggedLeftEnd ? sM->raggedStartStateProb
                                                          : sM->startStateProb); //Initialise forward matrix.

    //Backward matrix.
    DpMatrix *backwardDpMatrix = dpMatrix_construct(diagonalNumber, sM->stateNumber);

    int64_t tracedBackTo = 0;
    int64_t totalPosteriorCalculations = 0;
    while (1) { //Loop that moves through the matrix forward
        Diagonal diagonal = bandIterator_getNext(forwardBandIterator);

        //Forward calculation
        dpDiagonal_zeroValues(dpMatrix_createDiagonal(forwardDpMatrix, diagonal));
        diagonalCalculationForward(sM, diagonal_getXay(diagonal), forwardDpMatrix, sX, sY);

        bool atEnd = diagonal_getXay(diagonal) == diagonalNumber; //Condition true at the end of the matrix
        bool tracebackPoint = diagonal_getXay(diagonal) >= tracedBackTo + p->minDiagsBetweenTraceBack
                              && diagonal_getWidth(diagonal) <= p->diagonalExpansion * 2 +
                                                                1; //Condition true when we want to do an intermediate traceback.

        //Traceback
        if (atEnd || tracebackPoint) {
            //Initialise the last row (until now) of the backward matrix to represent an end point
            dpDiagonal_initialiseValues(dpMatrix_createDiagonal(backwardDpMatrix, diagonal), sM,
                                        (atEnd && alignmentHasRaggedRightEnd) ? sM->raggedEndStateProb
                                                                              : sM->endStateProb);
            if (diagonal_getXay(diagonal) >
                tracedBackTo + 1) { //This is a diagonal between the place we trace back to and where we trace back from
                DpDiagonal *j = dpMatrix_getDiagonal(forwardDpMatrix, diagonal_getXay(diagonal) - 1);
                assert(j != NULL);
                dpDiagonal_zeroValues(dpMatrix_createDiagonal(backwardDpMatrix, j->diagonal));
            }

            //Do walk back
            BandIterator *backwardBandIterator = bandIterator_clone(forwardBandIterator);
            Diagonal diagonal2 = bandIterator_getPrevious(backwardBandIterator);
            assert(diagonal_getXay(diagonal2) == diagonal_getXay(diagonal));
            int64_t tracedBackFrom = diagonal_getXay(diagonal) - (atEnd ? 0 : p->traceBackDiagonals + 1);
            double totalProbability = LOG_ZERO;
            int64_t totalPosteriorCalculationsThisTraceback = 0;
            while (diagonal_getXay(diagonal2) > tracedBackTo) {
                //Create the earlier diagonal
                if (diagonal_getXay(diagonal2) > tracedBackTo + 2) {
                    DpDiagonal *j = dpMatrix_getDiagonal(forwardDpMatrix, diagonal_getXay(diagonal2) - 2);
                    assert(j != NULL);
                    dpDiagonal_zeroValues(dpMatrix_createDiagonal(backwardDpMatrix, j->diagonal));
                }
                if (diagonal_getXay(diagonal2) > tracedBackTo + 1) {
                    diagonalCalculationBackward(sM, diagonal_getXay(diagonal2), backwardDpMatrix, sX, sY);
                }
                if (diagonal_getXay(diagonal2) <= tracedBackFrom) {
                    assert(dpMatrix_getDiagonal(forwardDpMatrix, diagonal_getXay(diagonal2)) != NULL);
                    assert(dpMatrix_getDiagonal(forwardDpMatrix, diagonal_getXay(diagonal2) - 1) != NULL);
                    assert(dpMatrix_getDiagonal(backwardDpMatrix, diagonal_getXay(diagonal2)) != NULL);
                    if (diagonal_getXay(diagonal2) != diagonalNumber) {
                        assert(dpMatrix_getDiagonal(backwardDpMatrix, diagonal_getXay(diagonal2) + 1) != NULL);
                    }
                    if (totalPosteriorCalculationsThisTraceback++ % 10 == 0) {
                        double newTotalProbability = diagonalCalculationTotalProbability(sM, diagonal_getXay(diagonal2),
                                                                                         forwardDpMatrix,
                                                                                         backwardDpMatrix, sX, sY);
                        if (totalPosteriorCalculationsThisTraceback != 1) {
                            if (totalProbability + 1.0 <= newTotalProbability || newTotalProbability + 1.0 <= newTotalProbability) {
                                char *logIdentifier = getLogIdentifier();
                                st_logInfo(" %s Error calculating newTotalProbability: totalProbability:%f, totalNewProbability:%f",
                                        logIdentifier, totalProbability, newTotalProbability);
                                free(logIdentifier);
                            }
                            assert(totalProbability + 1.0 > newTotalProbability);
                            assert(newTotalProbability + 1.0 > newTotalProbability);
                        }
                        totalProbability = newTotalProbability;
                    }

                    diagonalPosteriorProbFn(sM, diagonal_getXay(diagonal2), forwardDpMatrix, backwardDpMatrix, sX, sY,
                                            totalProbability, p, extraArgs);

                    if (diagonal_getXay(diagonal2) < tracedBackFrom || atEnd) {
                        dpMatrix_deleteDiagonal(forwardDpMatrix, diagonal_getXay(
                                diagonal2)); //Delete forward diagonal after last access in posterior calculation
                    }
                }
                if (diagonal_getXay(diagonal2) + 1 <= diagonalNumber) {
                    dpMatrix_deleteDiagonal(backwardDpMatrix, diagonal_getXay(diagonal2) +
                                                              1); //Delete backward diagonal after last access in backward calculation
                }
                diagonal2 = bandIterator_getPrevious(backwardBandIterator);
            }
            tracedBackTo = tracedBackFrom;
            bandIterator_destruct(backwardBandIterator);
            dpMatrix_deleteDiagonal(backwardDpMatrix, diagonal_getXay(diagonal2) + 1);
            dpMatrix_deleteDiagonal(forwardDpMatrix, diagonal_getXay(diagonal2));
            //Check memory state.
            assert(dpMatrix_getActiveDiagonalNumber(backwardDpMatrix) == 0);
            totalPosteriorCalculations += totalPosteriorCalculationsThisTraceback;
            if (!atEnd) {
                assert(dpMatrix_getActiveDiagonalNumber(forwardDpMatrix) == p->traceBackDiagonals + 2);
            }
        }

        if (atEnd) {
            break;
        }
    }
    assert(totalPosteriorCalculations == diagonalNumber);
    assert(tracedBackTo == diagonalNumber);
    assert(dpMatrix_getActiveDiagonalNumber(backwardDpMatrix) == 0);
    assert(dpMatrix_getActiveDiagonalNumber(forwardDpMatrix) == 0);
    //Cleanup
    dpMatrix_destruct(forwardDpMatrix);
    dpMatrix_destruct(backwardDpMatrix);
    bandIterator_destruct(forwardBandIterator);
    band_destruct(band);
}

/*
 * Computes for the forward log probability of aligning the two sequences
 */
double computeForwardProbability(SymbolString sX, SymbolString sY, stList *anchorPairs, PairwiseAlignmentParameters *p,
                                 StateMachine *sM,
                                 bool alignmentHasRaggedLeftEnd, bool alignmentHasRaggedRightEnd) {
    //Prerequisites
    assert(p->traceBackDiagonals >= 1);
    assert(p->diagonalExpansion >= 0);
    assert(p->diagonalExpansion % 2 == 0);
    assert(p->minDiagsBetweenTraceBack >= 2);
    assert(p->traceBackDiagonals + 1 < p->minDiagsBetweenTraceBack);

    int64_t diagonalNumber = sX.length + sY.length;
    if (diagonalNumber == 0) { //Deal with trivial case
        return LOG_ONE;
    }

    //Primitives for the forward matrix recursion
    Band *band = band_construct(anchorPairs, sX.length, sY.length, p->diagonalExpansion);
    BandIterator *forwardBandIterator = bandIterator_construct(band);
    DpMatrix *forwardDpMatrix = dpMatrix_construct(diagonalNumber, sM->stateNumber);
    dpDiagonal_initialiseValues(dpMatrix_createDiagonal(forwardDpMatrix, bandIterator_getNext(forwardBandIterator)), sM,
                                alignmentHasRaggedLeftEnd ? sM->raggedStartStateProb
                                                          : sM->startStateProb); //Initialise forward matrix.

    double totalLogProbability = LOG_ZERO;

    while (1) { //Loop that moves through the matrix forward
        Diagonal diagonal = bandIterator_getNext(forwardBandIterator);

        //Forward calculation
        dpDiagonal_zeroValues(dpMatrix_createDiagonal(forwardDpMatrix, diagonal));
        diagonalCalculationForward(sM, diagonal_getXay(diagonal), forwardDpMatrix, sX, sY);

        bool atEnd = diagonal_getXay(diagonal) == diagonalNumber; //Condition true at the end of the matrix
        if (atEnd) {
            //Backward matrix.
            DpMatrix *backwardDpMatrix = dpMatrix_construct(diagonalNumber, sM->stateNumber);
            dpDiagonal_initialiseValues(dpMatrix_createDiagonal(backwardDpMatrix, diagonal), sM,
                                        alignmentHasRaggedRightEnd ? sM->raggedEndStateProb : sM->endStateProb);
            totalLogProbability = diagonalCalculationTotalProbability(sM, diagonalNumber,
                                                                      forwardDpMatrix, backwardDpMatrix, sX, sY);
            dpMatrix_deleteDiagonal(backwardDpMatrix, diagonalNumber);
            dpMatrix_destruct(backwardDpMatrix);
            break;
        }
    }
    //Cleanup
    for (int64_t i = 0; i <= diagonalNumber; i++) {
        dpMatrix_deleteDiagonal(forwardDpMatrix, i);
    }
    dpMatrix_destruct(forwardDpMatrix);
    bandIterator_destruct(forwardBandIterator);
    band_destruct(band);

    return totalLogProbability;
}

///////////////////////////////////
///////////////////////////////////
//Split large gap functions
//
//Functions to split up alignment around gaps in the anchors that are too large.
///////////////////////////////////
///////////////////////////////////

static bool getSplitPointsP(int64_t *x1, int64_t *y1, int64_t x2, int64_t y2, int64_t x3, int64_t y3,
                            stList *splitPoints, int64_t splitMatrixBiggerThanThis, bool skipBlock) {
    /*
     * x2/y2 are the previous anchor point, x3/y3 are the next anchor point. Gaps greater than (x3-x2)*(y3-y2) are split up.
     */
    int64_t lX2 = x3 - x2;
    int64_t lY2 = y3 - y2;
    int64_t matrixSize = lX2 * lY2;
    if (matrixSize > splitMatrixBiggerThanThis) {
        st_logDebug("Split point found at x1: %" PRIi64 " x2: %" PRIi64 " y1: %" PRIi64 " y2: %" PRIi64 "\n", x2, x3,
                    y2, y3);
        int64_t maxSequenceLength = sqrt(splitMatrixBiggerThanThis);
        int64_t hX = lX2 / 2 > maxSequenceLength ? maxSequenceLength : lX2 / 2;
        int64_t hY = lY2 / 2 > maxSequenceLength ? maxSequenceLength : lY2 / 2;
        if (!skipBlock) {
            stList_append(splitPoints, stIntTuple_construct4(*x1, *y1, x2 + hX, y2 + hY));
        }
        *x1 = x3 - hX;
        *y1 = y3 - hY;
        return 1;
    }
    return 0;
}

stList *getSplitPoints(stList *anchorPairs, int64_t lX, int64_t lY, int64_t splitMatrixBiggerThanThis,
                       bool alignmentHasRaggedLeftEnd, bool alignmentHasRaggedRightEnd) {
    int64_t x1 = 0, y1 = 0, x2 = 0, y2 = 0;
    assert(lX >= 0);
    assert(lY >= 0);
    stList *splitPoints = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);
    for (int64_t i = 0; i < stList_length(anchorPairs); i++) {
        stIntTuple *anchorPair = stList_get(anchorPairs, i);
        int64_t x3 = stIntTuple_get(anchorPair, 0), y3 = stIntTuple_get(anchorPair, 1);
        getSplitPointsP(&x1, &y1, x2, y2, x3, y3, splitPoints, splitMatrixBiggerThanThis,
                        alignmentHasRaggedLeftEnd && i == 0);
        assert(x3 >= x2);
        assert(y3 >= y2);
        assert(x3 < lX);
        assert(y3 < lY);
        x2 = x3 + 1;
        y2 = y3 + 1;
    }
    if (!getSplitPointsP(&x1, &y1, x2, y2, lX, lY, splitPoints, splitMatrixBiggerThanThis,
                         alignmentHasRaggedLeftEnd && stList_length(anchorPairs) == 0) || !alignmentHasRaggedRightEnd) {
        stList_append(splitPoints, stIntTuple_construct4(x1, y1, lX, lY));
    }

    if (stList_length(splitPoints) > 1) {
        st_logDebug("For sequences of length %" PRIi64 " and %" PRIi64 " we got %" PRIi64 " splits\n", lX, lY,
                    stList_length(splitPoints));
    }
    return splitPoints;
}

static void convertAlignedPairs(stList *alignedPairs2, int64_t offsetX, int64_t offsetY) {
    /*
     * Convert the coordinates of the computed pairs.
     */
    for (int64_t k = 0; k < stList_length(alignedPairs2); k++) {
        stIntTuple *i = stList_get(alignedPairs2, k);
        assert(stIntTuple_length(i) == 3);
        stList_set(alignedPairs2, k,
                   stIntTuple_construct3(stIntTuple_get(i, 0), stIntTuple_get(i, 1) + offsetX,
                                         stIntTuple_get(i, 2) + offsetY));
        stIntTuple_destruct(i);
    }
}

void getPosteriorProbsWithBandingSplittingAlignmentsByLargeGaps(StateMachine *sM, stList *anchorPairs, SymbolString sX,
                                                                SymbolString sY,
                                                                PairwiseAlignmentParameters *p,
                                                                bool alignmentHasRaggedLeftEnd,
                                                                bool alignmentHasRaggedRightEnd,
                                                                void (*diagonalPosteriorProbFn)(StateMachine *, int64_t,
                                                                                                DpMatrix *, DpMatrix *,
                                                                                                const SymbolString,
                                                                                                const SymbolString,
                                                                                                double,
                                                                                                PairwiseAlignmentParameters *,
                                                                                                void *),
                                                                void (*coordinateCorrectionFn)(), void *extraArgs) {
    stList *splitPoints = getSplitPoints(anchorPairs, sX.length, sY.length, p->splitMatrixBiggerThanThis,
                                         alignmentHasRaggedLeftEnd, alignmentHasRaggedRightEnd);
    int64_t j = 0;
    //Now to the actual alignments
    for (int64_t i = 0; i < stList_length(splitPoints); i++) {
        stIntTuple *subRegion = stList_get(splitPoints, i);
        int64_t x1 = stIntTuple_get(subRegion, 0);
        int64_t y1 = stIntTuple_get(subRegion, 1);
        int64_t x2 = stIntTuple_get(subRegion, 2);
        int64_t y2 = stIntTuple_get(subRegion, 3);

        //Sub sequences
        SymbolString sX3 = symbolString_getSubString(sX, x1, x2 - x1);
        SymbolString sY3 = symbolString_getSubString(sY, y1, y2 - y1);

        //List of anchor pairs
        stList *subListOfAnchorPoints = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);
        while (j < stList_length(anchorPairs)) {
            stIntTuple *anchorPair = stList_get(anchorPairs, j);
            int64_t x = stIntTuple_get(anchorPair, 0);
            int64_t y = stIntTuple_get(anchorPair, 1);
            assert(x + y >= x1 + y1);
            if (x + y >= x2 + y2) {
                break;
            }
            assert(x >= x1 && x < x2);
            assert(y >= y1 && y < y2);
            stList_append(subListOfAnchorPoints, stIntTuple_construct3(x - x1, y - y1, stIntTuple_get(anchorPair, 2)));
            j++;
        }

        //Make the alignments
        getPosteriorProbsWithBanding(sM, subListOfAnchorPoints, sX3, sY3, p, (alignmentHasRaggedLeftEnd || i > 0),
                                     (alignmentHasRaggedRightEnd || i < stList_length(splitPoints) - 1),
                                     diagonalPosteriorProbFn, extraArgs);
        if (coordinateCorrectionFn != NULL) {
            coordinateCorrectionFn(x1, y1, extraArgs);
        }

        //Clean up
        stList_destruct(subListOfAnchorPoints);
        symbolString_destruct(sX3);
        symbolString_destruct(sY3);
    }
    assert(j == stList_length(anchorPairs));
    stList_destruct(splitPoints);
}

///////////////////////////////////
///////////////////////////////////
//Core public functions
///////////////////////////////////
///////////////////////////////////

PairwiseAlignmentParameters *pairwiseAlignmentBandingParameters_construct() {
    PairwiseAlignmentParameters *p = st_malloc(sizeof(PairwiseAlignmentParameters));
    p->threshold = 0.01;
    p->minDiagsBetweenTraceBack = 1000;
    p->traceBackDiagonals = 40;
    p->diagonalExpansion = 20;
    p->constraintDiagonalTrim = 14;
    p->splitMatrixBiggerThanThis = (int64_t) 3000 * 3000;
    p->alignAmbiguityCharacters = 0;
    p->gapGamma = 0.5;
    p->dynamicAnchorExpansion = 0;
    return p;
}

void pairwiseAlignmentBandingParameters_destruct(PairwiseAlignmentParameters *p) {
    free(p);
}

void *pairwiseAlignmentParameters_jsonParse(PairwiseAlignmentParameters *params, char *buf, size_t r) {
    // Setup parser
    jsmntok_t *tokens;
    char *js;
    int64_t tokenNumber = stJson_setupParser(buf, r, &tokens, &js);

    for (int64_t tokenIndex = 1; tokenIndex < tokenNumber; tokenIndex++) {
        jsmntok_t key = tokens[tokenIndex];
        char *keyString = stJson_token_tostr(js, &key);

        if (strcmp(keyString, "threshold") == 0) {
            params->threshold = stJson_parseFloat(js, tokens, ++tokenIndex);
        } else if (strcmp(keyString, "minDiagsBetweenTraceBack") == 0) {
            params->minDiagsBetweenTraceBack = stJson_parseInt(js, tokens, ++tokenIndex);
        } else if (strcmp(keyString, "traceBackDiagonals") == 0) {
            params->traceBackDiagonals = stJson_parseInt(js, tokens, ++tokenIndex);
        } else if (strcmp(keyString, "diagonalExpansion") == 0) {
            params->diagonalExpansion = stJson_parseInt(js, tokens, ++tokenIndex);
        } else if (strcmp(keyString, "constraintDiagonalTrim") == 0) {
            params->constraintDiagonalTrim = stJson_parseInt(js, tokens, ++tokenIndex);
        } else if (strcmp(keyString, "splitMatrixBiggerThanThis") == 0) {
            params->splitMatrixBiggerThanThis = stJson_parseInt(js, tokens, ++tokenIndex);
        } else if (strcmp(keyString, "alignAmbiguityCharacters") == 0) {
            params->alignAmbiguityCharacters = stJson_parseBool(js, tokens, ++tokenIndex);
        } else if (strcmp(keyString, "gapGamma") == 0) {
            params->gapGamma = stJson_parseFloat(js, tokens, ++tokenIndex);
        } else if (strcmp(keyString, "dynamicAnchorExpansion") == 0) {
            params->dynamicAnchorExpansion = stJson_parseBool(js, tokens, ++tokenIndex);
        } else {
            st_errAbort("ERROR: Unrecognised key in pairwise alignment parameters json: %s\n", keyString);
        }
    }

    // Cleanup
    free(js);
    free(tokens);
}

static void alignedPairCoordinateCorrectionFn(int64_t offsetX, int64_t offsetY, void *extraArgs) {
    stList *subListOfAlignedPairs = ((void **) extraArgs)[0];
    stList *alignedPairs = ((void **) extraArgs)[1];
    convertAlignedPairs(subListOfAlignedPairs, offsetX,
                        offsetY); //Shift back the aligned pairs to the appropriate coordinates
    while (stList_length(subListOfAlignedPairs) > 0) {
        stList_append(alignedPairs, stList_pop(subListOfAlignedPairs));
    }
}

static void pairCoordinateCorrectionFn(int64_t offsetX, int64_t offsetY, void *extraArgs) {
    for (int64_t i = 0; i < 6; i += 2) {
        stList *subListOfPairs = ((void **) extraArgs)[i];
        stList *pairs = ((void **) extraArgs)[i + 1];
        convertAlignedPairs(subListOfPairs, offsetX, offsetY); //Shift back the  pairs to the appropriate coordinates
        while (stList_length(subListOfPairs) > 0) {
            stList_append(pairs, stList_pop(subListOfPairs));
        }
    }
}

stList *getAlignedPairsUsingAnchors(StateMachine *sM, SymbolString sX, SymbolString sY, stList *anchorPairs,
                                    PairwiseAlignmentParameters *p,
                                    bool alignmentHasRaggedLeftEnd, bool alignmentHasRaggedRightEnd) {
    //This list of pairs to be returned. Not in any order, but points must be unique
    stList *subListOfAlignedPairs = stList_construct();
    stList *alignedPairs = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);
    void *extraArgs[2] = {subListOfAlignedPairs, alignedPairs};

    getPosteriorProbsWithBandingSplittingAlignmentsByLargeGaps(sM, anchorPairs, sX, sY, p,
                                                               alignmentHasRaggedLeftEnd, alignmentHasRaggedRightEnd,
                                                               diagonalCalculationPosteriorMatchProbs,
                                                               alignedPairCoordinateCorrectionFn, extraArgs);

    assert(stList_length(subListOfAlignedPairs) == 0);
    stList_destruct(subListOfAlignedPairs);

    return alignedPairs;
}

void getAlignedPairsWithIndelsUsingAnchors(StateMachine *sM, SymbolString sX, SymbolString sY, stList *anchorPairs,
                                           PairwiseAlignmentParameters *p, stList **alignedPairs, stList **gapXPairs,
                                           stList **gapYPairs,
                                           bool alignmentHasRaggedLeftEnd, bool alignmentHasRaggedRightEnd) {
    stList *subListOfAlignedPairs = stList_construct();
    *alignedPairs = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);

    stList *subListOfGapXPairs = stList_construct();
    *gapXPairs = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);

    stList *subListOfGapYPairs = stList_construct();
    *gapYPairs = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);

    void *extraArgs[6] = {subListOfAlignedPairs, *alignedPairs,
                          subListOfGapXPairs, *gapXPairs, subListOfGapYPairs, *gapYPairs};

    getPosteriorProbsWithBandingSplittingAlignmentsByLargeGaps(sM, anchorPairs, sX, sY, p,
                                                               alignmentHasRaggedLeftEnd, alignmentHasRaggedRightEnd,
                                                               diagonalCalculationPosteriorProbs,
                                                               pairCoordinateCorrectionFn, extraArgs);

    assert(stList_length(subListOfAlignedPairs) == 0);
    stList_destruct(subListOfAlignedPairs);
    assert(stList_length(subListOfGapXPairs) == 0);
    stList_destruct(subListOfGapXPairs);
    assert(stList_length(subListOfGapYPairs) == 0);
    stList_destruct(subListOfGapYPairs);
}

stList *getAlignedPairs(StateMachine *sM, SymbolString sX, SymbolString sY, PairwiseAlignmentParameters *p,
                        bool alignmentHasRaggedLeftEnd,
                        bool alignmentHasRaggedRightEnd) {
    stList *anchorPairs = stList_construct();
    stList *alignedPairs = getAlignedPairsUsingAnchors(sM, sX, sY, anchorPairs, p, alignmentHasRaggedLeftEnd,
                                                       alignmentHasRaggedRightEnd);
    stList_destruct(anchorPairs);
    return alignedPairs;
}

void getAlignedPairsWithIndels(StateMachine *sM, SymbolString sX, SymbolString sY, PairwiseAlignmentParameters *p,
                               stList **alignedPairs, stList **gapXPairs, stList **gapYPairs,
                               bool alignmentHasRaggedLeftEnd, bool alignmentHasRaggedRightEnd) {
    stList *anchorPairs = stList_construct();
    getAlignedPairsWithIndelsUsingAnchors(sM, sX, sY, anchorPairs, p,
                                          alignedPairs, gapXPairs, gapYPairs,
                                          alignmentHasRaggedLeftEnd, alignmentHasRaggedRightEnd);
    stList_destruct(anchorPairs);
}

void getExpectationsUsingAnchors(StateMachine *sM, Hmm *hmmExpectations, SymbolString sX, SymbolString sY,
                                 stList *anchorPairs,
                                 PairwiseAlignmentParameters *p, bool alignmentHasRaggedLeftEnd,
                                 bool alignmentHasRaggedRightEnd) {
    getPosteriorProbsWithBandingSplittingAlignmentsByLargeGaps(sM, anchorPairs, sX, sY, p,
                                                               alignmentHasRaggedLeftEnd, alignmentHasRaggedRightEnd,
                                                               diagonalCalculationExpectations, NULL,
                                                               hmmExpectations);
}

void getExpectations(StateMachine *sM, Hmm *hmmExpectations, SymbolString sX, SymbolString sY,
                     PairwiseAlignmentParameters *p,
                     bool alignmentHasRaggedLeftEnd, bool alignmentHasRaggedRightEnd) {
    stList *anchorPairs = stList_construct();
    getExpectationsUsingAnchors(sM, hmmExpectations, sX, sY, anchorPairs, p, alignmentHasRaggedLeftEnd,
                                alignmentHasRaggedRightEnd);
    stList_destruct(anchorPairs);
}

/*
 * Functions for adjusting weights to account for probability of alignment to a gap.
 */

int64_t *getIndelProbabilities(stList *alignedPairs, int64_t seqLength, bool xIfTrueElseY) {
    int64_t *indelProbs = st_malloc(seqLength * sizeof(int64_t));
    for (int64_t i = 0; i < seqLength; i++) {
        indelProbs[i] = PAIR_ALIGNMENT_PROB_1;
    }
    for (int64_t i = 0; i < stList_length(alignedPairs); i++) {
        stIntTuple *j = stList_get(alignedPairs, i);
        indelProbs[stIntTuple_get(j, xIfTrueElseY ? 1 : 2)] -= stIntTuple_get(j, 0);
    }
    for (int64_t i = 0; i < seqLength; i++) {
        if (indelProbs[i] < 0) {
            indelProbs[i] = 0;
        }
    }
    return indelProbs;
}

stList *reweightAlignedPairs(stList *alignedPairs,
                             int64_t *indelProbsX, int64_t *indelProbsY, double gapGamma) {
    stList *reweightedAlignedPairs = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);
    for (int64_t i = 0; i < stList_length(alignedPairs); i++) {
        stIntTuple *aPair = stList_get(alignedPairs, i);
        int64_t x = stIntTuple_get(aPair, 1);
        int64_t y = stIntTuple_get(aPair, 2);
        int64_t updatedWeight = stIntTuple_get(aPair, 0) - gapGamma * (indelProbsX[x] + indelProbsY[y]);
        stList_append(reweightedAlignedPairs, stIntTuple_construct3(updatedWeight, x, y));
    }
    stList_destruct(alignedPairs);
    return reweightedAlignedPairs;
}

stList *reweightAlignedPairs2(stList *alignedPairs, int64_t seqLengthX, int64_t seqLengthY, double gapGamma) {
    if (gapGamma <= 0.0) {
        return alignedPairs;
    }
    int64_t *indelProbsX = getIndelProbabilities(alignedPairs, seqLengthX, 1);
    int64_t *indelProbsY = getIndelProbabilities(alignedPairs, seqLengthY, 0);
    alignedPairs = reweightAlignedPairs(alignedPairs, indelProbsX, indelProbsY, gapGamma);
    free(indelProbsX);
    free(indelProbsY);
    return alignedPairs;
}

int64_t getNumberOfMatchingAlignedPairs(SymbolString subSeqX, SymbolString subSeqY, stList *alignedPairs) {
    int64_t matches = 0;
    for (int64_t i = 0; i < stList_length(alignedPairs); i++) {
        stIntTuple *aPair = stList_get(alignedPairs, i);
        int64_t x = stIntTuple_get(aPair, 1), y = stIntTuple_get(aPair, 2);
        matches += subSeqX.alphabet->symbolsEqual(subSeqX.sequence[x], subSeqY.sequence[y]) ? 1 : 0;
    }
    return matches;
}

double scoreByIdentity(SymbolString subSeqX, SymbolString subSeqY, stList *alignedPairs) {
    int64_t matches = getNumberOfMatchingAlignedPairs(subSeqX, subSeqY, alignedPairs);
    return 100.0 * ((subSeqX.length + subSeqY.length) == 0 ? 0 : (2.0 * matches) / (subSeqX.length + subSeqY.length));
}

double scoreByIdentityIgnoringGaps(SymbolString subSeqX, SymbolString subSeqY, stList *alignedPairs) {
    int64_t matches = getNumberOfMatchingAlignedPairs(subSeqX, subSeqY, alignedPairs);
    return 100.0 * matches / (double) stList_length(alignedPairs);
}

static double totalScore(stList *alignedPairs) {
    double score = 0.0;
    for (int64_t i = 0; i < stList_length(alignedPairs); i++) {
        stIntTuple *aPair = stList_get(alignedPairs, i);
        score += stIntTuple_get(aPair, 0);
    }
    return score;
}

double scoreByPosteriorProbability(int64_t lX, int64_t lY, stList *alignedPairs) {
    return 100.0 * ((lX + lY) == 0 ? 0 : (2.0 * totalScore(alignedPairs)) / ((lX + lY) * PAIR_ALIGNMENT_PROB_1));
}

double scoreByPosteriorProbabilityIgnoringGaps(stList *alignedPairs) {
    return 100.0 * totalScore(alignedPairs) / ((double) stList_length(alignedPairs) * PAIR_ALIGNMENT_PROB_1);
}

/*
 * Functions for pairwise alignment creation.
 */

static int64_t *getCumulativeGapProbs(stList *gapPairs, int64_t seqLength, bool seqXNotSeqY) {
    int64_t *gapCumulativeProbs = st_calloc(seqLength, sizeof(int64_t));

    // Work out the per-position gap probability
    for (int64_t i = 0; i < stList_length(gapPairs); i++) {
        stIntTuple *gapPair = stList_get(gapPairs, i);
        assert(stIntTuple_get(gapPair, seqXNotSeqY ? 1 : 2) >= 0);
        assert(stIntTuple_get(gapPair, seqXNotSeqY ? 1 : 2) < seqLength);
        gapCumulativeProbs[stIntTuple_get(gapPair, seqXNotSeqY ? 1 : 2)] += stIntTuple_get(gapPair, 0);
    }

    // Make cumulative
    for (int64_t i = 1; i < seqLength; i++) {
        gapCumulativeProbs[i] += gapCumulativeProbs[i - 1];
    }

    return gapCumulativeProbs;
}

static int64_t getIndelProb(int64_t *gapCumulativeProbs, int64_t start, int64_t length) {
    assert(start >= 0);
    assert(length >= 0);
    return length == 0 ? 0 : (gapCumulativeProbs[start + length - 1] - (start > 0 ? gapCumulativeProbs[start - 1] : 0));
}

stList *getMaximalExpectedAccuracyPairwiseAlignment(stList *alignedPairs,
                                                    stList *gapXPairs, stList *gapYPairs,
                                                    int64_t seqXLength, int64_t seqYLength, double *alignmentScore,
                                                    PairwiseAlignmentParameters *p) {

    int64_t totalPairs = stList_length(alignedPairs); // Total number of aligned pairs
    stList_sort(alignedPairs, cmpAlignedPairsByCoordinates);

    double *scores = st_calloc(totalPairs + 1, sizeof(double)); // MEA alignment score for each aligned pair
    int64_t *backPointers = st_calloc(totalPairs + 1, sizeof(int64_t)); // Trace back pointers
    bool *isHighScore = st_calloc(totalPairs + 1,
                                  sizeof(bool)); // Records if the score for a given aligned pair at index i is larger than any
    // score for an aligned pair at any index less than i

    // Calculate gap array cumulative probs
    int64_t *gapYCumulativeProbs = getCumulativeGapProbs(gapYPairs, seqYLength, FALSE);
    int64_t *gapXCumulativeProbs = getCumulativeGapProbs(gapXPairs, seqXLength, TRUE);

    // Iterate through the aligned pairs in order of increasing sequence coordinate

    double maxScore = 0; // Max score seen so far

    for (int64_t i = 0; i < totalPairs + 1; i++) {

        int64_t matchProb, x, y;

        if (i == totalPairs) { // Add final aligned pair at the end of the sequences to trace back the final alignment
            matchProb = 0;
            x = seqXLength;
            y = seqYLength;
        } else {
            stIntTuple *aPair = stList_get(alignedPairs, i);
            matchProb = stIntTuple_get(aPair, 0);
            x = stIntTuple_get(aPair, 1);
            y = stIntTuple_get(aPair, 2);
        }

        // The MEA alignment score of the pair with no preceding alignment pair
        double score = matchProb +
                       (getIndelProb(gapXCumulativeProbs, 0, x) + getIndelProb(gapYCumulativeProbs, 0, y)) *
                       p->gapGamma;
        int64_t backPointer = -1;

        // Walk back through previous aligned pairs
        for (int64_t j = i - 1; j >= 0; j--) {
            stIntTuple *pPair = stList_get(alignedPairs, j);
            int64_t x2 = stIntTuple_get(pPair, 1), y2 = stIntTuple_get(pPair, 2);

            // If the previous pair, pPair, and aPair can form an alignment
            if (x2 < x && y2 < y) {

                // Calc score of MEA alignment including pPair
                int64_t s = matchProb + scores[j] +
                            (getIndelProb(gapXCumulativeProbs, x2 + 1, x - x2 - 1) +
                             getIndelProb(gapYCumulativeProbs, y2 + 1, y - y2 - 1)) * p->gapGamma;

                // If score s is highest keep it
                if (s > score) {
                    score = s;
                    backPointer = j;
                }

                // If the score of pPair is a high score then can not increase score by exploring further back
                // pointers
                if (isHighScore[j]) {
                    break;
                }
            }
        }

        // Store the best alignment for aPair
        backPointers[i] = backPointer;
        scores[i] = score;

        // If the score of the alignment ending at aPair is higher than any we've seen to date
        double s = score + ((x < seqXLength ? getIndelProb(gapXCumulativeProbs, x + 1, seqXLength - x - 1) : 0) +
                            (y < seqYLength ? getIndelProb(gapYCumulativeProbs, y + 1, seqYLength - y - 1) : 0)) *
                           p->gapGamma;
        //fprintf(stderr, "scores %"PRId64"x, %"PRId64"y, %"PRId64"s, %"PRId64"S\n", x, y, score, maxScore);
        if (s >= maxScore) {
            maxScore = s; // Record the max score
            isHighScore[i] = 1; // Record the fact that the score represents a max seen so far.
        }
    }

    // Trace back to build the MEA alignment in reverse
    stList *filteredAlignment = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);
    int64_t i = backPointers[totalPairs];
    while (i >= 0) {
        stIntTuple *aPair = stList_get(alignedPairs, i);
        stList_append(filteredAlignment, stIntTuple_construct3(stIntTuple_get(aPair, 0),
                                                               stIntTuple_get(aPair, 1), stIntTuple_get(aPair, 2)));
        i = backPointers[i];
            }
    stList_reverse(filteredAlignment); // Flip the order

    // Cleanup
    free(scores);
    free(backPointers);
    free(isHighScore);
    free(gapXCumulativeProbs);
    free(gapYCumulativeProbs);

    *alignmentScore = maxScore;
    return filteredAlignment;
}

stList *leftShiftAlignment(stList *alignedPairs, SymbolString seqX, SymbolString seqY) {

    stList *leftShiftedAlignedPairs = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);

    int64_t x = seqX.length, y = seqY.length;
    for (int64_t i = stList_length(alignedPairs) - 1; i >= 0; i--) {
        stIntTuple *alignedPair = stList_get(alignedPairs, i);
        int64_t x2 = stIntTuple_get(alignedPair, 1), y2 = stIntTuple_get(alignedPair, 2);

        while ((x - x2 > 1 || y - y2 > 1) && seqX.alphabet->symbolsEqual(seqX.sequence[x - 1], seqY.sequence[y -
                                                                                                             1])) { // Insert in seqX or seqY and shift possible
            stList_append(leftShiftedAlignedPairs, stIntTuple_construct3(stIntTuple_get(alignedPair, 0), x - 1, y -
                                                                                                                1)); // Hacks the score by borrowing from the current aligned pair being considered
            x--;
            y--;

            if (x2 == x || y2 == y) { // We've shifted over an existing aligned pair
                break;
            }
        }
        if (x2 < x && y2 < y) {
            stList_append(leftShiftedAlignedPairs, stIntTuple_construct3(stIntTuple_get(alignedPair, 0), x2, y2));
            x = x2;
            y = y2;
        }
    }

    // Deal with boundary at beginning of alignment
    while ((x > 0 & y > 0) && seqX.alphabet->symbolsEqual(seqX.sequence[x - 1], seqY.sequence[y - 1])) {
        int64_t score = stList_length(alignedPairs) > 0 ? stIntTuple_get(stList_get(alignedPairs, 0), 0) : 1;
        stList_append(leftShiftedAlignedPairs, stIntTuple_construct3(score, x - 1, y - 1));
        x--;
        y--;
    }

    // Reverse, because built backwards
    stList_reverse(leftShiftedAlignedPairs);

    return leftShiftedAlignedPairs;
}

/*
 * Convenience function that aligns two sequences return a left-shift MEA alignment
 */
stList *
getShiftedMEAAlignment(SymbolString seqX, SymbolString seqY, stList *anchorAlignment, PairwiseAlignmentParameters *p,
                       StateMachine *sM,
                       bool alignmentHasRaggedLeftEnd, bool alignmentHasRaggedRightEnd, double *alignmentScore) {
    // Generate the posterior alignment probabilities
    stList *alignedPairs, *gapXPairs, *gapYPairs;
    getAlignedPairsWithIndelsUsingAnchors(sM, seqX, seqY, anchorAlignment,
                                          p, &alignedPairs, &gapXPairs, &gapYPairs,
                                          alignmentHasRaggedLeftEnd, alignmentHasRaggedRightEnd);

    // Get the MEA alignment
    stList *alignment = getMaximalExpectedAccuracyPairwiseAlignment(alignedPairs, gapXPairs, gapYPairs,
                                                                    seqX.length, seqY.length,
                                                                    alignmentScore, p);

    // Left shift the alignment
    stList *leftShiftedAlignment = leftShiftAlignment(alignment, seqX, seqY);

    // Cleanup
    stList_destruct(gapXPairs);
    stList_destruct(gapYPairs);
    stList_destruct(alignedPairs);
    stList_destruct(alignment);

    return leftShiftedAlignment;
}

stList *filterPairwiseAlignmentToMakePairsOrdered(stList *alignedPairs, SymbolString seqX, SymbolString seqY,
                                                  PairwiseAlignmentParameters *p) {
    stList *phonyGapPairs = stList_construct();
    double alignmentScore;
    stList *filteredAlignedPairs = getMaximalExpectedAccuracyPairwiseAlignment(alignedPairs, phonyGapPairs,
                                                                               phonyGapPairs,
                                                                               seqX.length, seqY.length,
                                                                               &alignmentScore, p);
    stList_destruct(phonyGapPairs);
    return filteredAlignedPairs;
}

/*
 * Code to create quick and dirty alignment anchors.
 */

#define KMER_SIZE 20 // that's an alphabet of 4^20 ~= 1x10^12 for DNA/RNA

uint64_t kmerKey(const void *k) {
    uint64_t hash = 0; //5381;
    Symbol *cA = (Symbol *) k;
    for (int64_t i = 0; i < KMER_SIZE; i++) {
        hash = cA[i] + (hash << 6) + (hash << 16) - hash;
    }
    return hash;
}

int kmerEqualKey(const void *key1, const void *key2) {
    Symbol *cA = (Symbol *) key1, *cA2 = (Symbol *) key2;
    for (int64_t i = 0; i < KMER_SIZE; i++) {
        if (cA[i] != cA2[i]) {
            return 0;
        }
    }
    return 1;
}

stHash *getKmers(SymbolString seq, uint64_t *l) {
    stHash *kmerOccurrences = stHash_construct3(kmerKey, kmerEqualKey, NULL, NULL);
    for (int64_t i = 0; i < seq.length - KMER_SIZE + 1; i++) {
        l[i] = i;
        stHash_insert(kmerOccurrences, &(seq.sequence[i]),
                      &(l[i])); // For speed and simplicity this only allows for one entry per k-mer
        // for KMER_SIZE >= 20, this should not be a big deal as collisions should be rare
    }
    return kmerOccurrences;
}

typedef struct _chainPair {
    uint64_t x, y, score;
    int64_t backpointer;
    bool highScore;
} ChainPair;

stList *getKmerAlignmentAnchors(SymbolString seqX, SymbolString seqY, uint64_t anchorExpansion) {
    if (KMER_SIZE > seqX.length || KMER_SIZE > seqY.length) { // Won't work if KMER_SIZE larger than sequences lengths
        return stList_construct();
    }

    // Get the kmers in seqX
    uint64_t l[seqX.length - KMER_SIZE + 1];
    stHash *kmerOccurrences = getKmers(seqX, l);

    // Get the kmers shared in both seqX and seqY
    ChainPair cPs[seqY.length - KMER_SIZE + 1]; // Array of chain pairs, representing shared k-mers
    uint64_t i = 0, maxScore = 0;
    int64_t maxPair = -1;
    for (int64_t y = 0; y < seqY.length - KMER_SIZE + 1; y++) {
        int64_t *x = stHash_search(kmerOccurrences, &(seqY.sequence[y]));
        if (x != NULL) {
            cPs[i].x = *x;
            cPs[i].y = y;
            cPs[i].score = 1; // score of chain pair if no prior chain-pair in chain
            cPs[i].backpointer = -1; // back-pointer to previous chain pair (-1 indicates no prior pair)

            // Find highest scoring subchain to extend
            for (int64_t j = i - 1; j >= 0; j--) { // walk through the previous pairs in decreasing y coordinate
                if (cPs[j].x <
                    cPs[i].x) { // if previous aligned k-mer pair has smaller x coordinate then can be chained
                    if (cPs[j].score + 1 > cPs[i].score) { // chain score is better than seen so far
                        cPs[i].score = cPs[j].score + 1; // update chain score
                        cPs[i].backpointer = j; // update backpointer
                    }
                    if (cPs[j].highScore) { // if this is a high score then can not do better by looking back further, so break
                        break;
                    }
                }
            }

            // Update if is equal to highest score seen so far
            if (cPs[i].score >= maxScore) {
                cPs[i].highScore = 1;
                maxScore = cPs[i].score;
                maxPair = i;
            } else {
                cPs[i].highScore = 0;
            }

            i++;
        }
    }

    // Traceback
    stList *anchorPairs = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);
    while (maxPair != -1) {
        int64_t x = cPs[maxPair].x + KMER_SIZE / 2;
        int64_t y = cPs[maxPair].y + KMER_SIZE / 2;
        stList_append(anchorPairs,
                      stIntTuple_construct3(x, y, anchorExpansion));
        maxPair = cPs[maxPair].backpointer;
    }
    stList_reverse(anchorPairs); // Put in ascending coordinate order

    // Cleanup
    stHash_destruct(kmerOccurrences);

    return anchorPairs;
}