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
#ifndef TYPHON_DISTRIBUTE_ELEMENT_CLOUD_H
#define TYPHON_DISTRIBUTE_ELEMENT_CLOUD_H

#include <string>

#include "typhon.h"



namespace _TYPH_Internal {

struct Displacements;
class Layer_Info;

class Element_Cloud
{
public:
    static int constexpr    CGEL = 0;
    static int constexpr    CREG = 1;
    static int constexpr    CMAT = 2;
    static int constexpr  COWNPE = 3;
    static int constexpr COWNLEL = 4;
    static int constexpr   CCONS = 5;

    explicit Element_Cloud(int nd_per_el);

    // Index data array
    int       &operator()(int i, int j);
    int const &operator()(int i, int j) const;

    void Set_Dims(int glob_size, int nproc, int rank);

    void Set_Data(int nproc, int rank, Layer_Info &layer_info, MPI_Comm *comm,
            int nel_glob, Displacements &sdispls, Displacements &rdispls);

    void Redistribute(int rank, Displacements &sdispls, Displacements &rdispls,
            Layer_Info &layer_info, MPI_Comm *comm);

    //void dump(std::string filename) const;

public:
    int const len;          // Number of properties in data array

    int      size =  0;
    int   my_size =  0;
    int my_offset = -1;

    int *data = nullptr;
};

} // namespace _TYPH_Internal



#endif // TYPHON_DISTRIBUTE_ELEMENT_CLOUD_H
