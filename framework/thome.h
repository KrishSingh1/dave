///////////////////////////////////////////////////////////////////////////////////////
/// \file thome.h
/// \brief Parse .csv files containing gridcell-level Thome values, used for
/// temperature acclimation of vcmax25.
///
/// \author Drew Holzworth
/// $Date: 2022-11-15 11:14:08 +1100 (Tue, 15 Nov 2022) $
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_THOME_H
#define LPJ_GUESS_THOME_H

#include <vector>
#include "gutil.h"

/*
Class for reading .csv files containing gridcell-level Thome values.
*/
class THomeInput {
public:
    /*
    Constructor.
    */
    THomeInput();

    /*
    Initialise the input class.

    @param file_path: Path to the input file.
    */
    void init(xtring file_path);

    /*
    Read the thome value for a specific gridcell.

    @param latitude: Latitude of the gridcell, in degrees.
    @param longitude: Longitude of the gridcell, in degrees.
    @return Long-term mean maximum temperature of warmest month of the gridcell.
    */
    double read_gridcell(double latitude, double longitude) const;
private:
    /*
    Represents a thome value in a particular gridcell.

    @param latitude: Latitude of the gridcell, in degrees.
    @param longitude: Longitude of the gridcell, in degrees.
    @param thome: Long-term mean maximum temperature of warmest month of the gridcell.
    */
    class Datum {
        public:
            /*
            Constructor.

            @param latitude: Latitude of the gridcell, in degrees.
            @param longitude: Longitude of the gridcell, in degrees.
            @param thome: Long-term mean maximum temperature of warmest month of the gridcell.
            */
            Datum(double lat, double lon, double t);

            // Get the latitude of the gridcell.
            double get_latitude() const;
            // Get the longitude of the gridcell.
            double get_longitude() const;
            // Get the thome value of the gridcell.
            double get_thome() const;
        private:
            // Latitude of the gridcell, in degrees.
            double latitude;
            // Longitude of the gridcell, in degrees.
            double longitude;
            // Long-term mean maximum temperature of warmest month of the gridcell.
            double thome;
    };

    /*
    True iff init() has been called.
    */
    bool is_initialised;
    /*
    The parsed Thome values.
    */
    std::vector<Datum> data;
};

#endif // LPJ_GUESS_THOME_H
