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
#include <iostream>
#include <fstream>
#include <unistd.h> // getopt

#include "typhon.h"
#include "types.h"
#include "serialise.h"

#ifdef TYPHON_YAMLCPP_SUPPORT
#include "yaml_output.h"
#endif



void
Usage(char *argv[])
{
    std::cerr << "Usage: " << argv[0] << " -o <format> <dump.typh>\n";
    std::cerr << "\n";
    std::cerr << "\t-o\toutput format (yaml)\n";
}



bool
Validate_Format(std::string format, char *argv[])
{
#ifdef TYPHON_YAMLCPP_SUPPORT
    if (format == "yaml") {
        return true;
    }
#endif

#if !(defined TYPHON_YAMLCPP_SUPPORT) // && ...
    std::cerr << "No output formats compiled in, please rebuild\n";
    return false;
#else
    if (format.length() > 0) {
        std::cerr << "Unrecognised output format '" << format << "'\n";
    } else {
        Usage(argv);
    }
    return false;
#endif
}



int
main(int argc, char *argv[])
{
    using namespace _TYPH_Internal;

    std::string format;

    int opt;
    while ((opt = getopt(argc, argv, "o:")) != -1) {
        switch (opt) {
            case 'o':
                format = std::string(optarg);
                break;

            default:
                Usage(argv);
                return EXIT_FAILURE;
        }
    }

    // Should be exactly 1 non-optional argument
    if (argc - optind != 1) {
        Usage(argv);
        return EXIT_FAILURE;
    }

    // Check that a valid output format was specified
    if (!Validate_Format(format, argv)) {
        return EXIT_FAILURE;
    }

    std::string filename(argv[optind]);
    std::ifstream ifs(filename.c_str(), std::ios::binary);
    if (!ifs.is_open()) {
        std::cerr << "Couldn't open " << filename << "\n";
        return EXIT_FAILURE;
    }

#ifdef TYPHON_YAMLCPP_SUPPORT
    //YAML::Node root;
    YAML::Emitter emitter;
    if (format == "yaml") {
        emitter << YAML::BeginSeq;
    }
#endif

    bool seen_proc_header = false;
    bool success = true;
    int marker;
    while (Get_Marker(ifs, &marker) == TYPH_SUCCESS) {

        switch (marker) {
        case PROC_HEADER_SERIALISATION_MARKER:
            {
                int rank, nproc;
                Read_Proc_Header(ifs, &rank, &nproc);

#ifdef TYPHON_YAMLCPP_SUPPORT
                if (format == "yaml") {
                    if (seen_proc_header) {
                        emitter << YAML::EndSeq;
                        emitter << YAML::EndMap;
                    }

                    emitter << YAML::BeginMap;
                    emitter << YAML::Key << "rank"  << YAML::Value << rank;
                    emitter << YAML::Key << "nproc" << YAML::Value << nproc;

                    emitter << YAML::Key << "data"
                            << YAML::Value << YAML::BeginSeq;
                }
#endif

                seen_proc_header = true;
            }
            break;

        case Key_Set::SERIALISATION_MARKER:
            {
                Key_Set *key_set = new Key_Set();
                Deserialise(ifs, key_set);
                if (!seen_proc_header) break;

#ifdef TYPHON_YAMLCPP_SUPPORT
                if (format == "yaml") emitter << *key_set;
#endif
            }
            break;

        case Phase::SERIALISATION_MARKER:
            {
                Phase *phase = new Phase();
                Deserialise(ifs, phase);
                if (!seen_proc_header) break;

#ifdef TYPHON_YAMLCPP_SUPPORT
                if (format == "yaml") emitter << *phase;
#endif
            }
            break;

        case Quant::SERIALISATION_MARKER:
            {
                Quant *quant = new Quant();
                Deserialise(ifs, quant);
                if (!seen_proc_header) break;

#ifdef TYPHON_YAMLCPP_SUPPORT
                if (format == "yaml") emitter << *quant;
#endif
            }
            break;

        default:
            std::cerr << "Unrecognised file " << filename << "\n";
            success = false;
        }

        if (!success) break;
    }

    ifs.close();

#ifdef TYPHON_YAMLCPP_SUPPORT
    if (format == "yaml") {
        if (seen_proc_header) {
            emitter << YAML::EndSeq; // proc data
            emitter << YAML::EndMap; // proc
            emitter << YAML::EndSeq; // root seq
        }
        std::cout << emitter.c_str();
    }
#endif

    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
