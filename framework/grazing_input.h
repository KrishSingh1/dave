////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \file grazing_input.h
/// \brief Class for grazing input data (used mainly for grazing).
/// File format is <start_date> <end_date> <num_cows??>
/// \author Drew Holzworth and Juergen Knauer
/// $Date$
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef GRAZING_INPUT_H
#define GRAZING_INPUT_H

#include "guess.h"

class GrazingDate {
public:
    GrazingDate();
    GrazingDate(int year, int month, int day);
    ~GrazingDate();
    int get_year() const;
    int get_month() const;
    int get_day() const;
    bool equals(const GrazingDate &other) const;
    bool less_than(const GrazingDate &other) const;
    bool less_than_or_equals(const GrazingDate &other) const;
    bool greater_than(const GrazingDate &other) const;
    bool greater_than_or_equals(const GrazingDate &other) const;
private:
    int year, month, day;
};

class GrazingRecord {
public:
    GrazingRecord(GrazingDate start, GrazingDate end, double animal_units);
    ~GrazingRecord();
    double get_animal_units() const;

    /*
    Check if the specified date lies within the window defined by this grazing
    record.
    */
    bool within(const GrazingDate &date) const;
private:
    GrazingDate start_date;
    GrazingDate end_date;
    double animal_units;
};

class GrazingInput {
public:
    GrazingInput();
    ~GrazingInput();

    /*
    Open the specified input file. Return true if the file was successfully
    opened.
    */
    bool init(xtring path);
    bool is_initialised() const;
    bool load(/*const Coord &c*/);

    /*
    Get the number of animal units present in the gridcell on this date.
    Returns 0 if no animals present on this date.
    */
    double get_animal_units(const GrazingDate &date) const;
private:
    std::vector<GrazingRecord> data;
    xtring file_path;
    bool initialised;

};

#endif // GRAZING_INPUT_H
