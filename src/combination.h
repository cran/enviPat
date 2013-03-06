//
//  combination.h
//  CalcIsoStruct
//
//  Created by Christian Gerber on 11/29/12.
//  Copyright (c) 2012 Eawag. All rights reserved.
//

#ifndef CalcIsoStruct_combination_h
#define CalcIsoStruct_combination_h

#include "isotope.h"
#include "element.h"

typedef struct Compound Compound;
typedef struct Combination Combination;
typedef struct CompoundMulti CompoundMulti;
typedef struct CombinationMulti CombinationMulti;

struct Compound{
    unsigned int sum[MAX_ISO_ELEM];
    double mass;
    double abundance;
    int counter;
};

struct Combination {
    Compound compounds[MAX_COMPOUNDS];
    Element element;
    double max_abundance;
    double max_mass;
    int amount;
};

struct CompoundMulti{
    unsigned short counter[MAX_ELEMENTS];
    int sum[MAX_ISO_SIZE];
    double mass;
    double abundance;
    unsigned short indicator_iso;
};

struct CombinationMulti {
    CompoundMulti compounds[MAX_COMPOUNDS];
    double max_abundance;
    double max_mass;
    int amount;
};

int create_combination(Combination* combination,
                       Element *element,
                       int n,
                       double threshold);
int compoundmulti_sort_by_abundance_dec(const void *a,
                                        const void *b);
int compoundmulti_sort_by_abundance_inc(const void *a,
                                        const void *b);
int create_combination_2(double* m,
                         double* a,
                         int *cc,
                         double* max_a,
                         Element *elements,
                         int element_amount,
                         double threshold,
                         unsigned int* peak_amount,
                         int peak_limit);

int calc_combinations(Combination* combinations,
                         double threshold,
                         unsigned short element_amount,
                         double* mass,
                         double* a,
                         unsigned int* peak_amount,
                         unsigned int peak_limit
                      );
int clean_combinations(Combination* combinations,
                       double threshold,
                       unsigned short comb_amount);
void print_compoundmulti(CompoundMulti cm);

#endif
