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
#ifndef TYPHON_DISTRIBUTE_LAYER_INFO_H
#define TYPHON_DISTRIBUTE_LAYER_INFO_H

#include <string>



namespace _TYPH_Internal {

class Layer_Info
{
public:
    Layer_Info();
    ~Layer_Info();

    // Index connectivity data
    int       &operator()(int i, int j);
    int const &operator()(int i, int j) const;

    //void Dump(std::string filename) const;

    void Allocate_Connectivity(int n1, int n2);
    void Deallocate_Connectivity();

public:
    int lay = -1;               // Layer index (0 <= lay <= num layers)

    int nel = 0;                // # elements in layer
    int nnd = 0;                // # nodes in layer
    int nd_per_el = 0;          // # nodes per element

    int eloff = 0;              // element offset
    int ndoff = 0;              // node offset

    int *el_glob = nullptr;     // Global element indices
    int *nd_glob = nullptr;     // Global node indices

    int *conn_data = nullptr;   // Connectivity data
    int *conn_dims = nullptr;   // Dimensions of connectivity data
};

struct Node_Cloud;

void
Get_Layer_Info(
        int nel_local,
        int nel_glob,
        int nnd_glob,
        int nd_per_el,
        int nproc,
        int rank,
        MPI_Comm *comm,
        int num_layers,
        int const *conn_data,
        int const *partition,
        Node_Cloud &nc,
        Layer_Info *layer_info_arr);

} // namespace _TYPH_Internal



#endif // TYPHON_DISTRIBUTE_LAYER_INFO_H
