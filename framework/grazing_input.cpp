#include <stdio.h>
#include <stdlib.h>

#include "grazing_input.h"

const char *START_DATE_COLUMN = "startdate";
const char *END_DATE_COLUMN = "enddate";
const char *ANIMAL_UNITS_COLUMN = "animal_unit_ha";


GrazingDate::GrazingDate() {
    year = month = day = 0;
}

GrazingDate::GrazingDate(int year, int month, int day) {
    this->year = year;
    this->month = month;
    this->day = day;

    assert(month > 0 && month < 13);
    assert(day > 0 && day < 32); // fixme
}

GrazingDate::~GrazingDate() {

}

int GrazingDate::get_year() const {
    return year;
}

int GrazingDate::get_month() const {
    return month;
}

int GrazingDate::get_day() const {
    return day;
}


bool GrazingDate::equals(const GrazingDate &other) const {
    return year == other.year && month == other.month && day == other.day;
}

bool GrazingDate::less_than(const GrazingDate &other) const {
    if (year > other.get_year()) {
        return false;
    }
    if (year < other.get_year()) {
        return true;
    }
    // Same year.
    if (month > other.get_month()) {
        return false;
    }
    if (month < other.get_month()) {
        return true;
    }
    // Same month.
    return day < other.get_day();
}

bool GrazingDate::greater_than(const GrazingDate &other) const {
    return !(less_than(other) || equals(other));
}

bool GrazingDate::less_than_or_equals(const GrazingDate &other) const {
    return less_than(other) || equals(other);
}

bool GrazingDate::greater_than_or_equals(const GrazingDate &other) const {
    return !less_than(other);
}

GrazingRecord::GrazingRecord(GrazingDate start, GrazingDate end, double animal_units) {
    start_date = start;
    end_date = end;
    this->animal_units = animal_units;
}

GrazingRecord::~GrazingRecord() {

}

bool GrazingRecord::within(const GrazingDate &date) const {
    return start_date.less_than_or_equals(date) && end_date.greater_than_or_equals(date);
}

double GrazingRecord::get_animal_units() const {
    return animal_units;
}

GrazingInput::GrazingInput() {
    initialised = false;
}

GrazingInput::~GrazingInput() {

}

char *get_field(char *line, int num, const char *delim) {
    // size_t delim_size = strlen(delim);
    // char *delim_newline = (char *)calloc(strlen(delim) + 2, sizeof(char));
    // sprintf(delim_newline, "%s\n", delim); // TODO Drew: account for CR characters

    char* tok;
    for (tok = strtok(line, delim); tok && *tok; tok = strtok(NULL, delim)) {
        if (!--num) {
            // free(delim_newline);
            return tok;
        }
    }
    // free(delim_newline);
    return NULL;
}

GrazingDate parse_date(const char *str) {
    // for now, expect dd/MM/yyyy.
    const char *delim = "/";

    char *tmp0 = strdup(str);
    char *tmp1 = strdup(str);
    char *tmp2 = strdup(str);

    const char *day_str   = get_field(tmp0, 1, delim);
    const char *month_str = get_field(tmp1, 2, delim);
    const char *year_str  = get_field(tmp2, 3, delim);

    int day = atoi(day_str);
    int month = atoi(month_str);
    int year = atoi(year_str);

    free(tmp0);
    free(tmp1);
    free(tmp2);

    GrazingDate date(year, month, day);
    return date;
}

bool GrazingInput::init(xtring path) {
    file_path = path;

    FILE *file = fopen((const char *)file_path, "r");
    if (!file) {
        fail("Failed to open grazing input file: '%s'", (const char *)path);
    }
    char line[1024];
    int i = 0;
    while (fgets(line, 1024, file))
    {
        const char *delim = ",\n";
        char *tmp0 = strdup(line);
        char *tmp1 = strdup(line);
        char *tmp2 = strdup(line);

        const char *start_date_str = get_field(tmp0, 1, ",\n");
        const char *end_date_str = get_field(tmp1, 2, ",\n");
        const char *animal_units_str = get_field(tmp2, 3, ",\n");

        if (i == 0) {
            // Header row - verify column names. // TODO Drew: account for BOM
            assert(!strcmp(start_date_str, START_DATE_COLUMN));
            assert(!strcmp(end_date_str, END_DATE_COLUMN));
            assert(!strcmp(animal_units_str, ANIMAL_UNITS_COLUMN));
        } else {
            // Non-header row - parse values.
            GrazingDate start_date = parse_date(start_date_str);
            GrazingDate end_date = parse_date(end_date_str);
            double animal_units = atof(animal_units_str);

            GrazingRecord record(start_date, end_date, animal_units);
            data.push_back(record);
        }

        free(tmp0);
        free(tmp1);
        free(tmp2);

        i++;
    }

    fclose(file);

    initialised = true;;
    return true;
}

bool GrazingInput::is_initialised() const {
    return initialised;
}

bool GrazingInput::load(/*const Coord &gridcell*/) {
    // TBI - think about how this will vary with location.
    return true;
}

double GrazingInput::get_animal_units(const GrazingDate &date) const {
    for (const GrazingRecord &record : data) {
        if (record.within(date)) {
            return record.get_animal_units();
        }
    }
    return 0;
}
