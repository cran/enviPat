//
//  combination.c
//  CalcIsoStruct
//
//  Created by Christian Gerber on 11/29/12.
//  Copyright (c) 2012 Eawag. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "element.h"
#include "combination.h"
#include "peak.h"
#include "parse.h"
#include "preferences.h"
#include "n-tuple.h"

int create_combination(Combination* combination, Element *element, int mass_amount, double threshold)
{
    
    double abundance = 1.0;
    double mass = 0.0;
    int c = 0;
    int last_c = 1;
    unsigned int sum_form[MAX_ISO_ELEM];
    
    for (int i  = 0; i< element->amount; i++) {
        mass += (element->isotopes)->mass;
        abundance *= (element->isotopes)->abundance;
    }
    
    combination->amount = 0;
    combination->max_abundance = abundance;
    combination->max_mass = mass;
    combination->element = *element;

    combination->compounds->mass = mass;
    combination->compounds->abundance =  abundance;
    combination->compounds->counter = 0;
    combination->compounds->sum[0] = element->amount;

    last_c = 1;
    c++;
    
    for (int j = 1; j < element->iso_amount; j++) {
        
        for (int l = 0; l < last_c; l++) {
            
            mass = (combination->compounds + l)->mass;
            abundance = (combination->compounds + l)->abundance;
            memcpy(sum_form, (combination->compounds + l)->sum, MAX_ISO_ELEM * sizeof(unsigned int));
            
            for (int k = (combination->compounds + l)->counter; k < element->amount; k++) {
                
                sum_form[j]++;
                
                mass -= element->isotopes->mass;
                mass += (element->isotopes + j)->mass;

                abundance *= (element->isotopes + j)->abundance * (sum_form[0]);
                abundance /= (element->isotopes )->abundance * sum_form[j];
                
                sum_form[0]--;

                if ( abundance * pow(mass_amount,element->iso_amount) *10000 > threshold ) {
                    (combination->compounds + c)->mass = mass;
                    (combination->compounds + c)->abundance =  abundance;
                    (combination->compounds + c)->counter = k + 1;
                    memcpy((combination->compounds + c)->sum, sum_form, MAX_ISO_ELEM * sizeof(unsigned int));
                    c++;
                }
                
                if (c > MAX_COMPOUNDS) {
                    return 1;
                }
                
                if (combination ->max_abundance < abundance) {
                    combination ->max_abundance = abundance;
                }
                
                if (combination ->max_mass < mass) {
                    combination ->max_mass = mass;
                }
            }
        }
        last_c = c - 1;
    }
    
    combination->amount =  c;
   	return 0;
}


int compoundmulti_sort_by_abundance_dec(const void *a, const void *b)
{
    double y1 = ((const struct CompoundMulti*)a)->abundance;
    double y2 = ((const struct CompoundMulti*)b)->abundance;
    
    if (y1 < y2) {
        return 1;
    } else if (y1 > y2) {
        return -1;
    }
    return 0;
}

int compoundmulti_sort_by_abundance_inc(const void *a, const void *b)
{
    double y1 = ((const struct CompoundMulti*)a)->abundance;
    double y2 = ((const struct CompoundMulti*)b)->abundance;
    
    if (y1 > y2) {
        return 1;
    } else if (y1 < y2) {
        return -1;
    }
    return 0;
}

void create_isotope_list(Element *elements, int element_amount, Isotope2 *isotopes, int *iso_c){
    
    *iso_c = 0;
    for (int i = 0; i < element_amount; i++) {
        for (int j = 1; j < (elements + i)->iso_amount; j++) {
            (isotopes + *iso_c)->element_nr = i;
            (isotopes + *iso_c)->iso_e_nr = j;
            (isotopes + *iso_c)->amount = (elements + i)->amount;
            (isotopes + *iso_c)->abundance = (((elements + i)->isotopes + j)->abundance);
            (isotopes + *iso_c)->mass = (((elements + i)->isotopes + j)->mass);
            strcpy((isotopes + *iso_c)->symbol, ((elements + i)->isotopes + j)->symbol);
            strcpy((isotopes + *iso_c)->isotope, ((elements + i)->isotopes + j)->isotope);
            *(iso_c) += 1;
        }
    }
    qsort(isotopes, *iso_c, sizeof(Isotope2), isotope2_sort_by_n_abundance_dec);
    
    for (int k = 0; k < *iso_c; k++) {
        (isotopes + k)->iso_nr = k;
    }
}

void calc_monoisotopic(Element* elements, int element_amount, CompoundMulti *monoisotopic){
    double abundance = 1.0;
    double mass = 0.0;
    unsigned short pos = 0;
    for (int i  = 0; i< element_amount; i++) {
        for (int g = 0; g< (elements + i)->amount; g++) {
            mass += ((elements + i)->isotopes)->mass;
            abundance *= ((elements + i)->isotopes)->abundance;
            monoisotopic->sum[pos]++;
        }
        pos+= (elements + i)->iso_amount;
        monoisotopic->counter[i] = 0;
    }
    monoisotopic->mass = mass;
    monoisotopic->abundance = abundance;
    monoisotopic->indicator_iso = 0;
}


int create_combination_2(double* m,
                         double* a,
                         int *cc,
                         double* max_a,
                         Element *elements,
                         int element_amount,
                         double threshold,
                         unsigned int* peak_amount,
                         int peak_limit)
{

    Isotope2 *isotopes = (Isotope2*)malloc(MAX_ISO_SIZE * sizeof(Isotope2));
    CompoundMulti* monoisotopic = (CompoundMulti*)calloc(1,sizeof(CompoundMulti));

    CombinationMulti* A = (CombinationMulti*)malloc(sizeof(CombinationMulti));
    CombinationMulti* A2 = (CombinationMulti*)malloc(sizeof(CombinationMulti));
    CombinationMulti* C = (CombinationMulti*)malloc(sizeof(CombinationMulti));
    
    int iso_c = 0;
    create_isotope_list(elements, element_amount, isotopes, &iso_c);
    calc_monoisotopic(elements, element_amount, monoisotopic);
    
    A->amount = 0;
    A->max_mass = 0.0;
    A->max_abundance = 1.0;
    
    A2->amount = 0;
    A2->max_mass = 0.0;
    A2->max_abundance = 1.0;

    C->amount = 0;
    C->max_mass = 0.0;
    C->max_abundance = 1.0;
    
    *m = monoisotopic->mass;
    *a = monoisotopic->abundance;
    *max_a = monoisotopic->abundance;
    memcpy(cc, monoisotopic->sum, MAX_ISO_SIZE * sizeof(int));
    
    CompoundMulti* current_highest = NULL;
    CompoundMulti* current = NULL;

    current = monoisotopic;
    unsigned short iso_pos[MAX_ELEMENTS];
    for (unsigned short d = 0; d < element_amount; d++) {
        iso_pos[d] = 0;
        for (unsigned short b = 0; b < d; b++) {
            iso_pos[d] += (elements + b)->iso_amount;
        }
    }
    
    unsigned short iso_nr_max = 0;
    unsigned int c = 1;
    unsigned short h = 0;
    while (current->abundance != -1.0) {
        current_highest = current;
        C->amount = 0;
        iso_nr_max = 0;
        
        for (unsigned short j = current->indicator_iso; j < iso_c; j++) {
            h = (isotopes + j)->element_nr;
            
            if ( current->counter[h] < (elements + h)->amount ) {

                Isotope *isotope = (elements + h)->isotopes;
                unsigned short iso_e_nr = (isotopes + j)->iso_e_nr;
                C->compounds[C->amount] = *current;
                CompoundMulti *comp = &C->compounds[C->amount];

                comp->counter[h]++;
                comp->indicator_iso = j;
                comp->sum[iso_pos[h] + iso_e_nr]++;
                
                comp->mass -= isotope->mass;
                comp->mass += (isotope + iso_e_nr)->mass;
                
                comp->abundance *= ( isotope + iso_e_nr)->abundance * (comp->sum[iso_pos[h]]);
                comp->abundance /= ( isotope )->abundance * comp->sum[iso_pos[h] + iso_e_nr];
                
                comp->sum[iso_pos[h]]--;
                
                if (current_highest->abundance < comp->abundance) {
                    current_highest = comp;
                }
                
                if (comp->abundance >= current->abundance) {
                    iso_nr_max = j;
                }
            
                C->amount++;
            }
        }
        
        if (c > peak_limit) {
            *peak_amount = c;
            free(A);
            free(A2);
            free(C);
            free(isotopes);
            free(monoisotopic);
            return 1;
        }

        if(current_highest->abundance > current->abundance){
            if (*max_a < current_highest->abundance) {
                *max_a = current_highest->abundance;
            }
            
            for (int v = C->amount - 1; v >= 0 ; v--) {
                if ( C->compounds[v].abundance != current_highest->abundance
                    ) {
                    if(C->compounds[v].indicator_iso < iso_nr_max) {
                        A->compounds[A->amount] = C->compounds[v];
                        A->amount++;
                        
                        if (A->amount > MAX_COMPOUNDS) {
                            *peak_amount = c;
                            free(A);
                            free(A2);
                            free(C);
                            free(isotopes);
                            free(monoisotopic);
                            return 1;
                        }
                        
                    }else{
                        if((100/ *max_a) * C->compounds[v].abundance >= threshold) {
                            
                            A2->compounds[A2->amount] = C->compounds[v];
                            A2->amount++;
                            
                            if (A2->amount > MAX_COMPOUNDS) {
                                *peak_amount = c;
                                free(A);
                                free(A2);
                                free(C);
                                free(isotopes);
                                free(monoisotopic);
                                return 1;
                            }
                        }
                    }
                }else{
                    C->compounds[C->amount] = C->compounds[v];
                    current_highest = &C->compounds[C->amount];
                }
            }

            current = current_highest;

            if ((100/ *max_a) * current_highest->abundance >= threshold) {
                    *(m + c) = current_highest->mass;
                    *(a + c) = current_highest->abundance;
                    memcpy((cc + c * MAX_ISO_SIZE), current_highest->sum, MAX_ISO_SIZE * sizeof(int));
                    c++;
            }
        }
        else{
            for (int v = 0; v < C->amount; v++) {
                if ( (100/ *max_a) * C->compounds[v].abundance >= threshold ) {
                    A2->compounds[A2->amount] = C->compounds[v];
                    A2->amount++;
                    
                    if (A2->amount > MAX_COMPOUNDS) {
                        *peak_amount = c;
                        free(A);
                        free(A2);
                        free(C);
                        free(isotopes);
                        free(monoisotopic);
                        return 1;
                    }
                }
            }
        
            if (A->amount > 0) {
                CompoundMulti *a_c = &A->compounds[A->amount - 1];
                current = a_c;
                
                if ((100/ *max_a) * a_c->abundance >= threshold) {
                    *(m + c) = a_c->mass;
                    *(a + c) = a_c->abundance;
                    memcpy((cc + c * MAX_ISO_SIZE), a_c->sum, MAX_ISO_SIZE * sizeof(int));
                    c++;
                }
                A->amount--;
            }else if(A->amount == 0 && A2->amount > 0){
                CompoundMulti *a2 = &A2->compounds[A2->amount - 1];
                current = a2;
                
                if ( (100/ *max_a) * a2->abundance >= threshold ) {
                    *(m + c) = a2->mass;
                    *(a + c) = a2->abundance;
                    memcpy((cc + c * MAX_ISO_SIZE), a2->sum, MAX_ISO_SIZE * sizeof(int));
                    c++;
                }
                A2->amount--;
                
            }else {
                current->abundance = -1.0;
            }
        }
    }
    *peak_amount =  c;
    
    free(A);
    free(A2);
    free(C);
    free(isotopes);
    free(monoisotopic);
   	return 0;
}


int calc_combinations(Combination* combinations,
                         double threshold,
                         unsigned short element_amount,
                         double* m,
                         double* a,
                         unsigned int* peak_amount,
                         unsigned int peak_limit){
    

    unsigned int ar[element_amount];
    const size_t len = element_amount;
    size_t sizes[element_amount];

    double max_abundance = 1.0;
    double max_mass = 0.0;
    
    for (int i = 0; i < element_amount; i++) {
        max_abundance *= (combinations + i)->max_abundance;
        max_mass += (combinations + i)->max_mass;
    }
    
    clean_combinations(combinations, threshold * max_abundance / 100 , element_amount);
    
    double mass;
    double abundance;
    for (int l = 0; l < element_amount; l++) {
        sizes[l] = (combinations + l)->amount;
        ar[l] = 0;
    }

    int v = 0;
    do {
        mass = 0.0;
        abundance = 1.0;
        
        for (int b = 0; b < element_amount; b++) {
            if (ar[b] <=  (combinations + b)->amount) {
                mass += ((combinations + b)->compounds + ar[b])->mass;
                abundance *= ((combinations + b)->compounds + ar[b])->abundance;
                
            }else{
                return 1;
            }
        }
        if (
            (100/max_abundance)* abundance >= threshold
            
        ) {
            *(m + v) = mass;
            *(a + v) = (100/max_abundance) * abundance;
            v++;
        }
        
        if (v > peak_limit) {
            break;
        }
    } while (MBnext_n_tuple(ar, len, sizes));
    *(peak_amount) = v;
    return 0;
}


int clean_combinations( Combination* combinations, double threshold, unsigned short comb_amount){

    double clean_abundance = 1.0;
    double clean_abundance_other = 1.0;
    for (unsigned short b = 0; b < comb_amount; b++) {
        clean_abundance_other = 1.0;
        for (int d = 0; d < comb_amount; d++) {
            if (d != b) {
                clean_abundance_other *= (combinations + d)->max_abundance;
            }
        }
        
        for (int c = (combinations + b)->amount - 1; c >= 0; c--) {
            clean_abundance = 1.0;
            clean_abundance *= ((combinations + b)->compounds + c)->abundance * clean_abundance_other;

            if ( clean_abundance < threshold) {
                if (c == (combinations + b)->amount - 1) {
                    (combinations + b)->amount--;
                    
                }else{
                    *((combinations + b)->compounds + c) = *((combinations + b)->compounds + (combinations + b)->amount  - 1);
                    (combinations + b)->amount--;
                }
            }
        }
    }
    
    return 0;
}


