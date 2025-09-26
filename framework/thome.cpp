#include <errno.h>
#include <stdio.h>
#include <stdlib.h>

#include "guessmath.h"
#include "shell.h"
#include "thome.h"

// Expected name of the latitude column in the input file.
const char *COL_LAT = "latitude";

// Expected name of the longitude column in the input file.
const char *COL_LON = "longitude";

// Expected name of the thome column in the input file.
const char *COL_THOME = "thome";

THomeInput::THomeInput() {
    is_initialised = false;
}

void THomeInput::init(xtring file_path) {
    data.clear();
    FILE *file = fopen((const char *)file_path, "r");
    if (!file) {
        fail("Failed to open %s for input: %s", (const char *)file_path, strerror(errno));
    }

    // Parse the file - expect CSV format with latitude, longitude, thome columns.
    // First row should be column headers.
    const size_t BUFSIZE = 1024;
    char buffer[BUFSIZE];
    if (!fgets(buffer, sizeof(buffer), file)) {
        fclose(file);
        fail("Failed to read header line from %s", (const char *)file_path);
    }

    // Parse data rows
    while (fgets(buffer, sizeof(buffer), file)) {
        char *token;
        double values[3];
        int i = 0;

        // Get first token (latitude)
        token = strtok(buffer, ",");
        while (token && i < 3) {
            char *endptr;
            values[i] = strtod(token, &endptr);
            if (endptr == token || *endptr != '\0' && *endptr != '\n' && *endptr != '\r') {
                fclose(file);
                fail("Invalid numeric value in %s", (const char *)file_path);
            }
            token = strtok(NULL, ",");
            i++;
        }

        if (i != 3) {
            fclose(file);
            fail("Invalid number of columns in %s", (const char *)file_path);
        }

        // Create and store new Datum object
        data.push_back(Datum(values[0], values[1], values[2]));
    }

    // Check for read errors
    if (ferror(file)) {
        fclose(file);
        fail("Error reading from %s: %s", (const char *)file_path, strerror(errno));
    }

    int err = fclose(file);
    if (err != 0) {
        fail("Failed to close %s: %s", (const char *)file_path, strerror(errno));
    }
    is_initialised = true;
}

double THomeInput::read_gridcell(double latitude, double longitude) const {
    if (!is_initialised) {
        fail("Attempted to call THomeInput::read_gridcell() before initialising");
    }
    for (size_t i = 0; i < data.size(); i++) {
        if (equal(latitude, data[i].get_latitude()) &&
            equal(longitude, data[i].get_longitude())) {
            return data[i].get_thome();
        }
    }
    fail("Failed to find thome value for gridcell at (%f, %f)", latitude, longitude);
    return -1; // Cannot happen.
}

THomeInput::Datum::Datum(double lat, double lon, double t) {
    latitude = lat;
    longitude = lon;
    thome = t;
}

double THomeInput::Datum::get_thome() const {
    return thome;
}

double THomeInput::Datum::get_latitude() const {
    return latitude;
}

double THomeInput::Datum::get_longitude() const {
    return longitude;
}
