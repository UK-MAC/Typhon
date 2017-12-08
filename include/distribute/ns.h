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
#ifndef TYPHON_DISTRIBUTE_NS_H
#define TYPHON_DISTRIBUTE_NS_H



namespace _TYPH_Internal {

// TODO(timrlaw):
//      What exactly does NS stand for?
//      This whole class needs sorting out. It works but it's very tangled and
//      unclear.
//
class NS
{
public:
    static int constexpr NS_W = 4;
    static int constexpr NS_N = NS_W + 1;

    void Init(int size);
    void Cleanup();
    void Resize(int new_top);
    bool New(int ipd, int indx);

private:
    int *p = nullptr;
    int top = 0;
    int size = 0;
    int *d = nullptr;
    int p_size = 0;
};

} // namespace _TYPH_Internal



#endif // TYPHON_DISTRIBUTE_NS_H
