//
//  element.c
//  CalcIsoStruct
//
//  Created by Christian Gerber on 11/29/12.
//  Copyright (c) 2012 Eawag. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "element.h"
#include "isotope.h"

int set_element(Element* element, Isotope* isotopes, char* name, int amount, int iso_amount)
{
    element->amount = amount;
    element->iso_amount = iso_amount;
    strcpy(element->name, name);
    memcpy(element->isotopes, isotopes, iso_amount*sizeof(Isotope));
    qsort(element->isotopes, iso_amount, sizeof(Isotope), isotope_sort_by_abundance);
	return 0;
}
