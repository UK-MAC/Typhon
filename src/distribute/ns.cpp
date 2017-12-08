/* @HEADER@
 * Crown Copyright 2018 AWE.
 *
 * This file is part of Typhon.
 *
 * Typhon is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 * 
 * Typhon is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * Typhon. If not, see http://www.gnu.org/licenses/.
 * @HEADER@ */
#include "typhon.h"
#include "distribute/ns.h"



namespace _TYPH_Internal {

auto constexpr IXns = &Index_2D<NS::NS_N>;

void
NS::Init(int _size)
{
    // Initialise to store references from 0 to size-1
    p = new int[_size];
    std::fill(p, p + _size, -1);
    p_size = _size;

    top  = -1;
    size = 0;
}



void
NS::Cleanup()
{
    if (p != nullptr) {
        delete[] p;
        p = nullptr;
    }

    if (d != nullptr) {
        delete[] d;
        d = nullptr;
    }

    top  = -1;
    size = 0;
}



void
NS::Resize(int new_top)
{
    int constexpr NS_INCR = 10000;

    // Resize d storage size to be at least new_top (actually new_top+NS_INCR)
    int const new_size = new_top + NS_INCR;
    int *tmp = new int[NS_N * new_size];

    // Copy old entries over
    if (top > -1) {
        for (int i = 0; i < top + 1; i++) {
            for (int j = 0; j < NS_N; j++) {
                tmp[IXns(j, i)] = d[IXns(j, i)];
            }
        }

        delete[] d;
    }

    // Nullify the rest
    for (int i = top + 1; i < new_size; i++) {
        for (int j = 0; j < NS_N; j++) {
            tmp[IXns(j, i)] = -1;
        }
    }

    d = tmp;
    size = new_size;
}



bool
NS::New(int ipd, int indx)
{
    int nproc;
    TYPH_Get_Size(&nproc);

    // Input PE number:  0 <= ipd  < nproc
    // Reference number: 0 <= indx < SIZE(ns_p) - 1 (-1??)
    assert(0 <= ipd  && ipd  < nproc);
    assert(0 <= indx && indx < p_size);

    // Returns false if reference has been sent to PE ip before, and true
    // otherwise. In the latter case, the transfer top PE ip is recorded, so
    // that further calls with the same arguments will return false.
    // TODO(timrlaw)
    //ip=ipd+1_ink ! Actually store ip+1, so can use zero as the null value

    bool ret = true;
    if (p[indx] == -1) {
        // There have been no sends for this reference, get an entry in ns_d to
        // store transfer. Increase size if necessary
        if (top + 1 >= size) {
            Resize(top + 2);
        }

        top++;
        p[indx] = top;

        // Record transfer to ip is first entry
        d[IXns(0, top)] = ipd;

    } else {
        // This is where the sorted list of PEs starts for this reference
        int ii = p[indx];

        while (true) {
            // Search for entry for ip
            int k;
            for (k = 0; k < NS_W; k++) {
                if ((d[IXns(k, ii)] == -1) || (ipd <= d[IXns(k, ii)])) {
                    break;
                }
            }

            // Check if have not gone off the end of ns_d(:,ii)
            if (k < NS_W) {
                if (ipd == d[IXns(k, ii)]) {
                    // Entry found, return false
                    ret = false;

                } else {
                    // Entry for ip not present, will return true
                    // Add entry for ip and shift remaining entries forward
                    int itmp = ipd;
                    while (true) {
                        int kk;
                        for (kk = k; kk < NS_W; kk++) {
                            int itmp2 = d[IXns(kk, ii)];
                            d[IXns(kk, ii)] = itmp;
                            if (itmp2 == -1) break;
                            itmp = itmp2;
                        }

                        if (kk < NS_W) return ret; // Found null entry, done

                        // Reached end of ns_d(:,ii)
                        if (d[IXns(NS_W, ii)] == -1) {

                            // No further rows linked on.  Link to new row, with
                            // single entry ip.  Increase size of ns_d if
                            // necessary.
                            if (top + 1 >= size) {
                                Resize(top + 2);
                            }

                            top++;
                            d[IXns(0, top)] = itmp;
                            d[IXns(NS_W, ii)] = top;
                            break; // done
                        }

                        // Follow link to start of new row
                        ii = d[IXns(NS_W, ii)];
                        k = 0;
                    }
                }

                return ret; // done
            } // k < NS_W

            // Reached the end of ns_d(:,ii)
            if (d[IXns(NS_W, ii)] == -1) {

                // Will return true, no further rows linked on. Link to new row,
                // with single entry ip, increase size of ns_d if necessary.
                if (top + 1 >= size) {
                    Resize(top + 2);
                }

                top++;
                d[IXns(0, top)] = ipd;
                d[IXns(NS_W, ii)] = top;
                break;
            }

            // Follow link to start of new row
            ii = d[IXns(NS_W, ii)];

        } // while (true)
    } // p[indx] != -1

    return ret;
}

} // namespace _TYPH_Internal
