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
#ifndef TYPHON_DISTRIBUTE_NODE_CLOUD_H
#define TYPHON_DISTRIBUTE_NODE_CLOUD_H

#include <string>

#include "typhon.h"



namespace _TYPH_Internal {

class NS;
class Layer_Info;

struct Node_Cloud
{
public:
    static int constexpr    CGNOD = 0;  // Global node number
    static int constexpr  COWNLEL = 1;  // Local element number
    static int constexpr     CGEL = 2;  // Global element number
    static int constexpr   COWNPE = 3;  // Owner rank
    static int constexpr COWNLNOD = 4;  // Local node number

    static int constexpr LEN = 5;   // constant # properties

    ~Node_Cloud();

    void Set_Dims(int glob_size, int npes, int my_rank);

    void Set_Data(int nproc, int rank, MPI_Comm *comm,
            Layer_Info const &layer_info, int nnd_glob, int *&nd_to_el,
            int *nd_to_el_dims, NS *ns);

    //void dump(std::string filename, Error &err) const;

public:
    int         size =  0;
    int      my_size =  0;
    int    my_offset = -1;
    int my_data_size =  0;

    int *data = nullptr;
    int *indx = nullptr;
    int *data_ref = nullptr;
};

} // namespace _TYPH_Internal



#endif // TYPHON_DISTRIBUTE_NODE_CLOUD_H
