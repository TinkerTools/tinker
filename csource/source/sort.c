/* sort.f -- translated by f2c (version 20050501).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ######################################################### */
/*     ##                                                     ## */
/*     ##  subroutine sort  --  heapsort of an integer array  ## */
/*     ##                                                     ## */
/*     ######################################################### */


/*     "sort" takes an input list of integers and sorts it */
/*     into ascending order using the Heapsort algorithm */


/* Subroutine */ int sort_(integer *n, integer *list)
{
    static integer i__, j, k, index, lists;



/*     perform the heapsort of the input list */

    /* Parameter adjustments */
    --list;

    /* Function Body */
    k = *n / 2 + 1;
    index = *n;
    while(*n > 1) {
	if (k > 1) {
	    --k;
	    lists = list[k];
	} else {
	    lists = list[index];
	    list[index] = list[1];
	    --index;
	    if (index <= 1) {
		list[1] = lists;
		return 0;
	    }
	}
	i__ = k;
	j = k + k;
	while(j <= index) {
	    if (j < index) {
		if (list[j] < list[j + 1]) {
		    ++j;
		}
	    }
	    if (lists < list[j]) {
		list[i__] = list[j];
		i__ = j;
		j += j;
	    } else {
		j = index + 1;
	    }
	}
	list[i__] = lists;
    }
    return 0;
} /* sort_ */



/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine sort2  --  heapsort of real array with keys  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "sort2" takes an input list of reals and sorts it */
/*     into ascending order using the Heapsort algorithm; */
/*     it also returns a key into the original ordering */


/* Subroutine */ int sort2_(integer *n, doublereal *list, integer *key)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j, k, keys, index;
    static doublereal lists;



/*     initialize index into the original ordering */

    /* Parameter adjustments */
    --key;
    --list;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	key[i__] = i__;
    }

/*     perform the heapsort of the input list */

    k = *n / 2 + 1;
    index = *n;
    while(*n > 1) {
	if (k > 1) {
	    --k;
	    lists = list[k];
	    keys = key[k];
	} else {
	    lists = list[index];
	    keys = key[index];
	    list[index] = list[1];
	    key[index] = key[1];
	    --index;
	    if (index <= 1) {
		list[1] = lists;
		key[1] = keys;
		return 0;
	    }
	}
	i__ = k;
	j = k + k;
	while(j <= index) {
	    if (j < index) {
		if (list[j] < list[j + 1]) {
		    ++j;
		}
	    }
	    if (lists < list[j]) {
		list[i__] = list[j];
		key[i__] = key[j];
		i__ = j;
		j += j;
	    } else {
		j = index + 1;
	    }
	}
	list[i__] = lists;
	key[i__] = keys;
    }
    return 0;
} /* sort2_ */



/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine sort3  --  heapsort of integer array with keys  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "sort3" takes an input list of integers and sorts it */
/*     into ascending order using the Heapsort algorithm; */
/*     it also returns a key into the original ordering */


/* Subroutine */ int sort3_(integer *n, integer *list, integer *key)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j, k, keys, index, lists;



/*     initialize index into the original ordering */

    /* Parameter adjustments */
    --key;
    --list;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	key[i__] = i__;
    }

/*     perform the heapsort of the input list */

    k = *n / 2 + 1;
    index = *n;
    while(*n > 1) {
	if (k > 1) {
	    --k;
	    lists = list[k];
	    keys = key[k];
	} else {
	    lists = list[index];
	    keys = key[index];
	    list[index] = list[1];
	    key[index] = key[1];
	    --index;
	    if (index <= 1) {
		list[1] = lists;
		key[1] = keys;
		return 0;
	    }
	}
	i__ = k;
	j = k + k;
	while(j <= index) {
	    if (j < index) {
		if (list[j] < list[j + 1]) {
		    ++j;
		}
	    }
	    if (lists < list[j]) {
		list[i__] = list[j];
		key[i__] = key[j];
		i__ = j;
		j += j;
	    } else {
		j = index + 1;
	    }
	}
	list[i__] = lists;
	key[i__] = keys;
    }
    return 0;
} /* sort3_ */



/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine sort4  --  heapsort of integer absolute values  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "sort4" takes an input list of integers and sorts it into */
/*     ascending absolute value using the Heapsort algorithm */


/* Subroutine */ int sort4_(integer *n, integer *list)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k, index, lists;



/*     perform the heapsort of the input list */

    /* Parameter adjustments */
    --list;

    /* Function Body */
    k = *n / 2 + 1;
    index = *n;
    while(*n > 1) {
	if (k > 1) {
	    --k;
	    lists = list[k];
	} else {
	    lists = list[index];
	    list[index] = list[1];
	    --index;
	    if (index <= 1) {
		list[1] = lists;
		return 0;
	    }
	}
	i__ = k;
	j = k + k;
	while(j <= index) {
	    if (j < index) {
		if ((i__1 = list[j], abs(i__1)) < (i__2 = list[j + 1], abs(
			i__2))) {
		    ++j;
		}
	    }
	    if (abs(lists) < (i__1 = list[j], abs(i__1))) {
		list[i__] = list[j];
		i__ = j;
		j += j;
	    } else {
		j = index + 1;
	    }
	}
	list[i__] = lists;
    }
    return 0;
} /* sort4_ */



/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine sort5  --  heapsort of integer array modulo m  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "sort5" takes an input list of integers and sorts it */
/*     into ascending order based on each value modulo "m" */


/* Subroutine */ int sort5_(integer *n, integer *list, integer *m)
{
    static integer i__, j, k, jmod, smod, j1mod, index, lists;



/*     perform the heapsort of the input list */

    /* Parameter adjustments */
    --list;

    /* Function Body */
    k = *n / 2 + 1;
    index = *n;
    while(*n > 1) {
	if (k > 1) {
	    --k;
	    lists = list[k];
	} else {
	    lists = list[index];
	    list[index] = list[1];
	    --index;
	    if (index <= 1) {
		list[1] = lists;
		return 0;
	    }
	}
	i__ = k;
	j = k + k;
	while(j <= index) {
	    if (j < index) {
		jmod = list[j] % *m;
		j1mod = list[j + 1] % *m;
		if (jmod < j1mod) {
		    ++j;
		} else if (jmod == j1mod && list[j] < list[j + 1]) {
		    ++j;
		}
	    }
	    smod = lists % *m;
	    jmod = list[j] % *m;
	    if (smod < jmod) {
		list[i__] = list[j];
		i__ = j;
		j += j;
	    } else if (smod == jmod && lists < list[j]) {
		list[i__] = list[j];
		i__ = j;
		j += j;
	    } else {
		j = index + 1;
	    }
	}
	list[i__] = lists;
    }
    return 0;
} /* sort5_ */



/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  subroutine sort6  --  heapsort of a text string array  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "sort6" takes an input list of character strings and sorts */
/*     it into alphabetical order using the Heapsort algorithm */


/* Subroutine */ int sort6_(integer *n, char *list, ftnlen list_len)
{
    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__, j, k, index;
    static char lists[256];


#define list_ref(a_0,a_1) &list[(a_1)*list_len + a_0]



/*     perform the heapsort of the input list */

    /* Parameter adjustments */
    list -= list_len;

    /* Function Body */
    k = *n / 2 + 1;
    index = *n;
    while(*n > 1) {
	if (k > 1) {
	    --k;
	    s_copy(lists, list_ref(0, k), (ftnlen)256, list_len);
	} else {
	    s_copy(lists, list_ref(0, index), (ftnlen)256, list_len);
	    s_copy(list_ref(0, index), list_ref(0, 1), list_len, list_len);
	    --index;
	    if (index <= 1) {
		s_copy(list_ref(0, 1), lists, list_len, (ftnlen)256);
		return 0;
	    }
	}
	i__ = k;
	j = k + k;
	while(j <= index) {
	    if (j < index) {
		if (s_cmp(list_ref(0, j), list_ref(0, j + 1), list_len, 
			list_len) < 0) {
		    ++j;
		}
	    }
	    if (s_cmp(lists, list_ref(0, j), (ftnlen)256, list_len) < 0) {
		s_copy(list_ref(0, i__), list_ref(0, j), list_len, list_len);
		i__ = j;
		j += j;
	    } else {
		j = index + 1;
	    }
	}
	s_copy(list_ref(0, i__), lists, list_len, (ftnlen)256);
    }
    return 0;
} /* sort6_ */

#undef list_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine sort7  --  heapsort of text strings with keys  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "sort7" takes an input list of character strings and sorts it */
/*     into alphabetical order using the Heapsort algorithm; it also */
/*     returns a key into the original ordering */


/* Subroutine */ int sort7_(integer *n, char *list, integer *key, ftnlen 
	list_len)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__, j, k, keys, index;
    static char lists[256];


#define list_ref(a_0,a_1) &list[(a_1)*list_len + a_0]



/*     initialize index into the original ordering */

    /* Parameter adjustments */
    --key;
    list -= list_len;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	key[i__] = i__;
    }

/*     perform the heapsort of the input list */

    k = *n / 2 + 1;
    index = *n;
    while(*n > 1) {
	if (k > 1) {
	    --k;
	    s_copy(lists, list_ref(0, k), (ftnlen)256, list_len);
	    keys = key[k];
	} else {
	    s_copy(lists, list_ref(0, index), (ftnlen)256, list_len);
	    keys = key[index];
	    s_copy(list_ref(0, index), list_ref(0, 1), list_len, list_len);
	    key[index] = key[1];
	    --index;
	    if (index <= 1) {
		s_copy(list_ref(0, 1), lists, list_len, (ftnlen)256);
		key[1] = keys;
		return 0;
	    }
	}
	i__ = k;
	j = k + k;
	while(j <= index) {
	    if (j < index) {
		if (s_cmp(list_ref(0, j), list_ref(0, j + 1), list_len, 
			list_len) < 0) {
		    ++j;
		}
	    }
	    if (s_cmp(lists, list_ref(0, j), (ftnlen)256, list_len) < 0) {
		s_copy(list_ref(0, i__), list_ref(0, j), list_len, list_len);
		key[i__] = key[j];
		i__ = j;
		j += j;
	    } else {
		j = index + 1;
	    }
	}
	s_copy(list_ref(0, i__), lists, list_len, (ftnlen)256);
	key[i__] = keys;
    }
    return 0;
} /* sort7_ */

#undef list_ref




/*     ######################################################### */
/*     ##                                                     ## */
/*     ##  subroutine sort8  --  heapsort to unique integers  ## */
/*     ##                                                     ## */
/*     ######################################################### */


/*     "sort8" takes an input list of integers and sorts it into */
/*     ascending order using the Heapsort algorithm, duplicate */
/*     values are removed from the final sorted list */


/* Subroutine */ int sort8_(integer *n, integer *list)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j, k, index, lists;



/*     perform the heapsort of the input list */

    /* Parameter adjustments */
    --list;

    /* Function Body */
    k = *n / 2 + 1;
    index = *n;
    while(*n > 1) {
	if (k > 1) {
	    --k;
	    lists = list[k];
	} else {
	    lists = list[index];
	    list[index] = list[1];
	    --index;
	    if (index <= 1) {
		list[1] = lists;

/*     remove duplicate values from final list */

		j = 1;
		i__1 = *n;
		for (i__ = 2; i__ <= i__1; ++i__) {
		    if (list[i__ - 1] != list[i__]) {
			++j;
			list[j] = list[i__];
		    }
		}
		if (j < *n) {
		    *n = j;
		}
		return 0;
	    }
	}
	i__ = k;
	j = k + k;
	while(j <= index) {
	    if (j < index) {
		if (list[j] < list[j + 1]) {
		    ++j;
		}
	    }
	    if (lists < list[j]) {
		list[i__] = list[j];
		i__ = j;
		j += j;
	    } else {
		j = index + 1;
	    }
	}
	list[i__] = lists;
    }
    return 0;
} /* sort8_ */



/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  subroutine sort9  --  heapsort to unique real values  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     "sort9" takes an input list of reals and sorts it into */
/*     ascending order using the Heapsort algorithm, duplicate */
/*     values are removed from the final sorted list */


/* Subroutine */ int sort9_(integer *n, doublereal *list)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j, k, index;
    static doublereal lists;



/*     perform the heapsort of the input list */

    /* Parameter adjustments */
    --list;

    /* Function Body */
    k = *n / 2 + 1;
    index = *n;
    while(*n > 1) {
	if (k > 1) {
	    --k;
	    lists = list[k];
	} else {
	    lists = list[index];
	    list[index] = list[1];
	    --index;
	    if (index <= 1) {
		list[1] = lists;

/*     remove duplicate values from final list */

		j = 1;
		i__1 = *n;
		for (i__ = 2; i__ <= i__1; ++i__) {
		    if (list[i__ - 1] != list[i__]) {
			++j;
			list[j] = list[i__];
		    }
		}
		if (j < *n) {
		    *n = j;
		}
		return 0;
	    }
	}
	i__ = k;
	j = k + k;
	while(j <= index) {
	    if (j < index) {
		if (list[j] < list[j + 1]) {
		    ++j;
		}
	    }
	    if (lists < list[j]) {
		list[i__] = list[j];
		i__ = j;
		j += j;
	    } else {
		j = index + 1;
	    }
	}
	list[i__] = lists;
    }
    return 0;
} /* sort9_ */



/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine sort10  --  heapsort to unique text strings  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "sort10" takes an input list of character strings and sorts */
/*     it into alphabetical order using the Heapsort algorithm, */
/*     duplicate values are removed from the final sorted list */


/* Subroutine */ int sort10_(integer *n, char *list, ftnlen list_len)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__, j, k, index;
    static char lists[256];


#define list_ref(a_0,a_1) &list[(a_1)*list_len + a_0]



/*     perform the heapsort of the input list */

    /* Parameter adjustments */
    list -= list_len;

    /* Function Body */
    k = *n / 2 + 1;
    index = *n;
    while(*n > 1) {
	if (k > 1) {
	    --k;
	    s_copy(lists, list_ref(0, k), (ftnlen)256, list_len);
	} else {
	    s_copy(lists, list_ref(0, index), (ftnlen)256, list_len);
	    s_copy(list_ref(0, index), list_ref(0, 1), list_len, list_len);
	    --index;
	    if (index <= 1) {
		s_copy(list_ref(0, 1), lists, list_len, (ftnlen)256);

/*     remove duplicate values from final list */

		j = 1;
		i__1 = *n;
		for (i__ = 2; i__ <= i__1; ++i__) {
		    if (s_cmp(list_ref(0, i__ - 1), list_ref(0, i__), 
			    list_len, list_len) != 0) {
			++j;
			s_copy(list_ref(0, j), list_ref(0, i__), list_len, 
				list_len);
		    }
		}
		if (j < *n) {
		    *n = j;
		}
		return 0;
	    }
	}
	i__ = k;
	j = k + k;
	while(j <= index) {
	    if (j < index) {
		if (s_cmp(list_ref(0, j), list_ref(0, j + 1), list_len, 
			list_len) < 0) {
		    ++j;
		}
	    }
	    if (s_cmp(lists, list_ref(0, j), (ftnlen)256, list_len) < 0) {
		s_copy(list_ref(0, i__), list_ref(0, j), list_len, list_len);
		i__ = j;
		j += j;
	    } else {
		j = index + 1;
	    }
	}
	s_copy(list_ref(0, i__), lists, list_len, (ftnlen)256);
    }
    return 0;
} /* sort10_ */

#undef list_ref


