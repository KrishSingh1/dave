///////////////////////////////////////////////////////////////////////////////////////
/// \file indata_test.cpp
/// \brief Compatibility tests for useful text input behavior ahead of its rewrite.
///
/// This Source Code Form is subject to the terms of the Mozilla Public
/// License, v. 2.0. If a copy of the MPL was not distributed with this
/// file, You can obtain one at http://mozilla.org/MPL/2.0/.
///////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"
#include <filesystem>
#include "indata.h"
#include "parameters.h"
#include "temporary_directory.h"

#include <array>
#include <fstream>
#include <stdexcept>
#include <string>

namespace {

using TextInput::TimeDataD;
using TextInput::fileformat;
using Row = std::array<double, 3>;

const int first_year = 2000;
const char* const column_names[] = {"A", "B", "C"};
const std::array<std::array<double, 2>, 3> coordinates = {{
    {{12.25, -33.75}}, {{-120.25, 45.25}}, {{150.25, -20.25}}
}};
const std::array<std::array<Row, 3>, 3> expected = {{
    {{{{0, -2.5, 0.001}}, {{4.25, 5, 600}}, {{7.75, -8, 9.125}}}},
    {{{{10.5, 11.25, 12}}, {{13, 0, -15.5}}, {{16.25, 17.5, 18}}}},
    {{{{20, -21.25, 22.5}}, {{23.75, 24, 25.25}}, {{26.5, 27.75, 0}}}}
}};

// Literal fixtures and independent expected values make shifted columns and
// cell/year mixups visible. No cache layout, search heuristic, or string-column
// naming convention is part of this contract.
const char* const global_static =
    "static A B C\n"
    "0 0 -2.5 1e-3\n";
const char* const global_yearly_rows =
    "2000 0 -2.5 1e-3\n"
    "2001 4.25 5 6E+2\n"
    "2002 7.75 -8 9.125\n";
const char* const local_static_rows =
    "12.25 -33.75 0 -2.5 1e-3\n"
    "-120.25 45.25 10.5 11.25 12\n"
    "150.25 -20.25 20 -21.25 22.5\n";
const char* const local_yearly_rows =
    "12.25 -33.75 2000 0 -2.5 1e-3\n"
    "12.25 -33.75 2001 4.25 5 6E+2\n"
    "12.25 -33.75 2002 7.75 -8 9.125\n"
    "-120.25 45.25 2000 10.5 11.25 12\n"
    "-120.25 45.25 2001 13 0 -15.5\n"
    "-120.25 45.25 2002 16.25 17.5 18\n"
    "150.25 -20.25 2000 20 -21.25 22.5\n"
    "150.25 -20.25 2001 23.75 24 25.25\n"
    "150.25 -20.25 2002 26.5 27.75 0\n";
const char* const local_yearly_headerless =
    "12.25 -33.75\n"
    "2000 0 -2.5 1e-3\n"
    "2001 4.25 5 6E+2\n"
    "2002 7.75 -8 9.125\n"
    "-120.25 45.25\n"
    "2000 10.5 11.25 12\n"
    "2001 13 0 -15.5\n"
    "2002 16.25 17.5 18\n"
    "150.25 -20.25\n"
    "2000 20 -21.25 22.5\n"
    "2001 23.75 24 25.25\n"
    "2002 26.5 27.75 0\n";

class InputFile {
public:
    explicit InputFile(const std::string& contents)
        : path((std::filesystem::path(directory.get_path()) / "input.txt").string()) {
        std::ofstream output(path.c_str(), std::ios::binary);
        output << contents;
        output.close();
        if (!output) throw std::runtime_error("Could not write text input fixture");
    }

    // Declare the file before its reader so input and sidecars are removed only
    // after the reader has released its handles, including on assertion failure.
    TemporaryDirectory directory;
    const std::string path;
};

class ScopedMapTextFile {
public:
    explicit ScopedMapTextFile(bool enabled) : previous(map_text_file) {
        map_text_file = enabled;
    }
    ~ScopedMapTextFile() { map_text_file = previous; }
    ScopedMapTextFile(const ScopedMapTextFile&) = delete;
    ScopedMapTextFile& operator=(const ScopedMapTextFile&) = delete;
private:
    bool previous;
};

enum class Access { Direct, Mapped, Memory };

Coord coordinate(std::size_t cell) {
    Coord result;
    result.id = static_cast<int>(cell);
    result.lon = coordinates[cell][0];
    result.lat = coordinates[cell][1];
    return result;
}

void add_coordinate(ListArray_id<Coord>& gridlist, std::size_t cell) {
    Coord& result = gridlist.createobj();
    result.id = static_cast<int>(cell);
    result.lon = coordinates[cell][0];
    result.lat = coordinates[cell][1];
}

bool open_local(TimeDataD& input, const InputFile& file, Access access) {
    if (access == Access::Direct) return input.Open(file.path.c_str());
    ScopedMapTextFile mapping(access == Access::Mapped);
    ListArray_id<Coord> gridlist;
    for (std::size_t cell = 0; cell < coordinates.size(); ++cell)
        add_coordinate(gridlist, cell);
    return input.Open(file.path.c_str(), gridlist);
}

void check_row(const TimeDataD& input, int year, const Row& values, bool named = true) {
    for (int column = 0; column < 3; ++column) {
        INFO("year=" << year << ", column=" << column);
        CHECK(input.Get(year, column) == Approx(values[column]));
        if (named)
            CHECK(input.Get(year, column_names[column]) == Approx(values[column]));
    }
}

void check_bulk(const TimeDataD& input, std::size_t cell, bool yearly) {
    const int years = yearly ? 3 : 1;
    std::array<double, 9> all;
    all.fill(-12345);
    input.Get(all.data());
    for (int year = 0; year < years; ++year) {
        Row row = {{-12345, -12345, -12345}};
        input.Get(first_year + year, row.data());
        for (int column = 0; column < 3; ++column) {
            INFO("cell=" << cell << ", year=" << year << ", column=" << column);
            CHECK(row[column] == Approx(expected[cell][year][column]));
            CHECK(all[3 * year + column] == Approx(expected[cell][year][column]));
            CHECK(row[column] == Approx(input.Get(first_year + year, column)));
        }
    }
}

void check_names(const TimeDataD& input) {
    std::array<std::array<char, TextInput::MAXNAMESIZE>, TextInput::MAXRECORDS> storage = {};
    std::array<char*, TextInput::MAXRECORDS> names;
    for (std::size_t i = 0; i < names.size(); ++i) names[i] = storage[i].data();
    REQUIRE(input.GetHeader(names.data()));
    for (int column = 0; column < 3; ++column) {
        REQUIRE(input.GetHeader(column) != nullptr);
        CHECK(std::string(input.GetHeader(column)) == column_names[column]);
        CHECK(std::string(names[column]) == column_names[column]);
    }
}

void check_local(const InputFile& file, fileformat format, Access access) {
    TimeDataD input;
    REQUIRE(open_local(input, file, access));
    const bool yearly = format == TextInput::LOCAL_YEARLY;
    CHECK(input.GetFormat() == format);
    CHECK(input.GetNCells() == 3);
    CHECK(input.GetnColumns() == 3);
    CHECK(input.GetnYears() == (yearly ? 3 : 1));
    if (yearly) CHECK(input.GetFirstyear() == first_year);
    check_names(input);
    for (std::size_t cell : {2, 0, 1, 2}) {
        INFO("cell=" << cell);
        REQUIRE(input.Load(coordinate(cell)));
        CHECK(input.isloaded());
        for (int year = 0; year < (yearly ? 3 : 1); ++year)
            check_row(input, first_year + year, expected[cell][year]);
        if (!yearly) {
            check_row(input, -1000, expected[cell][0]);
            check_row(input, 3000, expected[cell][0]);
        }
        // Memory-backed coordinate and bulk access have dedicated defect tests.
        if (access != Access::Memory) {
            CHECK(input.GetCoord().lon == Approx(coordinates[cell][0]));
            CHECK(input.GetCoord().lat == Approx(coordinates[cell][1]));
            check_bulk(input, cell, yearly);
        }
    }
}

} // namespace

TEST_CASE("TimeDataD reads headed global static data for every coordinate and year", "[timedata]") {
    InputFile file(global_static);
    TimeDataD input;
    REQUIRE(input.Open(file.path.c_str()));
    CHECK(input.GetFormat() == TextInput::GLOBAL_STATIC);
    CHECK(input.GetNCells() == 1);
    CHECK(input.GetnColumns() == 3);
    CHECK(input.GetnYears() == 1);
    CHECK(input.isloaded());
    check_names(input);
    for (std::size_t cell : {2, 0, 1}) {
        REQUIRE(input.Load(coordinate(cell)));
        for (int year : {-1000, 2000, 3000}) check_row(input, year, expected[0][0]);
    }
    check_bulk(input, 0, false);
}

TEST_CASE("TimeDataD reads global yearly data with and without a header", "[timedata]") {
    bool headed = true;
    SECTION("named columns") {}
    SECTION("explicit format and positional columns") { headed = false; }
    InputFile file((headed ? "Year A B C\n" : std::string()) + global_yearly_rows);
    TimeDataD input(headed ? TextInput::EMPTY : TextInput::GLOBAL_YEARLY);
    REQUIRE(input.Open(file.path.c_str()));
    CHECK(input.GetFormat() == TextInput::GLOBAL_YEARLY);
    CHECK(input.GetNCells() == 1);
    CHECK(input.GetnColumns() == 3);
    CHECK(input.GetnYears() == 3);
    CHECK(input.GetFirstyear() == first_year);
    CHECK(input.isloaded());
    if (headed) check_names(input);
    for (std::size_t cell : {2, 0, 1}) {
        REQUIRE(input.Load(coordinate(cell)));
        for (int year = 0; year < 3; ++year)
            check_row(input, first_year + year, expected[0][year], headed);
    }
    check_bulk(input, 0, true);
}

TEST_CASE("TimeDataD selects headed local static data in arbitrary cell order", "[timedata]") {
    InputFile file(std::string("Lon Lat A B C\n") + local_static_rows);
    Access access = Access::Direct;
    SECTION("direct file access") {}
    SECTION("mapped file access") { access = Access::Mapped; }
    SECTION("memory access") { access = Access::Memory; }
    check_local(file, TextInput::LOCAL_STATIC, access);
}

TEST_CASE("TimeDataD selects headed local yearly data in arbitrary cell order", "[timedata]") {
    InputFile file(std::string("Lon Lat Year A B C\n") + local_yearly_rows);
    Access access = Access::Direct;
    SECTION("direct file access") {}
    SECTION("mapped file access") { access = Access::Mapped; }
    SECTION("memory access") { access = Access::Memory; }
    check_local(file, TextInput::LOCAL_YEARLY, access);
}

TEST_CASE("TimeDataD reads headerless local yearly coordinate blocks", "[timedata]") {
    InputFile file(local_yearly_headerless);
    TimeDataD input(TextInput::LOCAL_YEARLY);
    REQUIRE(input.Open(file.path.c_str()));
    CHECK(input.GetFormat() == TextInput::LOCAL_YEARLY);
    CHECK(input.GetNCells() == 3);
    CHECK(input.GetnColumns() == 3);
    CHECK(input.GetnYears() == 3);
    CHECK(input.GetFirstyear() == first_year);
    for (std::size_t cell : {2, 0, 1, 2}) {
        REQUIRE(input.Load(coordinate(cell)));
        CHECK(input.isloaded());
        for (int year = 0; year < 3; ++year)
            check_row(input, first_year + year, expected[cell][year], false);
        check_bulk(input, cell, true);
    }
}

TEST_CASE("TimeDataD holds endpoint values outside the supplied yearly period", "[timedata]") {
    // This is the compatibility default, not a restriction on future policies.
    Access access = Access::Direct;
    bool global = false;
    SECTION("global") { global = true; }
    SECTION("local direct") {}
    SECTION("local mapped") { access = Access::Mapped; }
    SECTION("local memory") { access = Access::Memory; }
    InputFile file(global ? std::string("Year A B C\n") + global_yearly_rows
                          : std::string("Lon Lat Year A B C\n") + local_yearly_rows);
    TimeDataD input;
    REQUIRE(open_local(input, file, access));
    for (std::size_t cell : {2, 0, 1}) {
        REQUIRE(input.Load(coordinate(cell)));
        const std::size_t source = global ? 0 : cell;
        check_row(input, 1999, expected[source][0]);
        check_row(input, 1000, expected[source][0]);
        check_row(input, 2003, expected[source][2]);
        check_row(input, 3000, expected[source][2]);
    }
}

TEST_CASE("TimeDataD accesses named values independently of file column order", "[timedata]") {
    InputFile ordered("Year A B C\n2000 1 2 3\n2001 4 5 6\n");
    InputFile reordered("Year C A B\n2000 3 1 2\n2001 6 4 5\n");
    TimeDataD first, second;
    REQUIRE(first.Open(ordered.path.c_str()));
    REQUIRE(second.Open(reordered.path.c_str()));
    REQUIRE(second.GetHeader(0) != nullptr);
    CHECK(std::string(second.GetHeader(0)) == "C");
    CHECK(second.Get(2000, 0) == Approx(3));
    for (int year = 2000; year <= 2001; ++year)
        for (const char* name : column_names)
            CHECK(first.Get(year, name) == Approx(second.Get(year, name)));
}

TEST_CASE("TimeDataD distinguishes an absent column from a present zero", "[timedata]") {
    InputFile file(std::string("Lon Lat Year A B C\n") + local_yearly_rows);
    Access access = Access::Direct;
    SECTION("direct") {}
    SECTION("mapped") { access = Access::Mapped; }
    SECTION("memory") { access = Access::Memory; }
    TimeDataD input;
    REQUIRE(open_local(input, file, access));
    REQUIRE(input.Load(coordinate(0)));
    char present[] = "A", absent[] = "optional";
    CHECK(input.item_in_header(present));
    CHECK_FALSE(input.item_in_header(absent));
    CHECK(input.Get(2000, "A") == Approx(0));
    CHECK(input.Get(2000, "optional", true) == TextInput::NOTFOUND);
}

TEST_CASE("TimeDataD handles headed files containing a single cell or year", "[timedata]") {
    std::string contents;
    fileformat format = TextInput::EMPTY;
    int cells = 1, years = 1;
    SECTION("global single year") {
        contents = "Year A B C\n2000 1 2 3\n";
        format = TextInput::GLOBAL_YEARLY;
    }
    SECTION("local static single cell") {
        contents = "Lon Lat A B C\n12.25 -33.75 1 2 3\n";
        format = TextInput::LOCAL_STATIC;
    }
    SECTION("local yearly single cell") {
        contents = "Lon Lat Year A B C\n12.25 -33.75 2000 1 2 3\n12.25 -33.75 2001 4 5 6\n";
        format = TextInput::LOCAL_YEARLY;
        years = 2;
    }
    SECTION("local yearly single year") {
        contents = "Lon Lat Year A B C\n12.25 -33.75 2000 1 2 3\n-120.25 45.25 2000 4 5 6\n";
        format = TextInput::LOCAL_YEARLY;
        cells = 2;
    }
    InputFile file(contents);
    TimeDataD input;
    REQUIRE(input.Open(file.path.c_str()));
    CHECK(input.GetFormat() == format);
    CHECK(input.GetNCells() == cells);
    CHECK(input.GetnYears() == years);
    REQUIRE(input.Load(coordinate(0)));
    check_row(input, 2000, {{1, 2, 3}});
    if (years == 2) check_row(input, 2001, {{4, 5, 6}});
    if (cells == 2) {
        REQUIRE(input.Load(coordinate(1)));
        check_row(input, 2000, {{4, 5, 6}});
    }
}

TEST_CASE("TimeDataD sequential access returns each local cell once", "[timedata]") {
    bool yearly = true, headed = true;
    SECTION("headed yearly") {}
    SECTION("headed static") { yearly = false; }
    SECTION("headerless yearly") { headed = false; }
    InputFile file(!headed ? local_yearly_headerless
                          : yearly ? std::string("Lon Lat Year A B C\n") + local_yearly_rows
                                   : std::string("Lon Lat A B C\n") + local_static_rows);
    TimeDataD input(headed ? TextInput::EMPTY : TextInput::LOCAL_YEARLY);
    REQUIRE(input.Open(file.path.c_str()));
    for (std::size_t cell = 0; cell < 3; ++cell) {
        REQUIRE(input.LoadNext());
        CHECK(input.GetCoord().lon == Approx(coordinates[cell][0]));
        CHECK(input.GetCoord().lat == Approx(coordinates[cell][1]));
        for (int year = 0; year < (yearly ? 3 : 1); ++year)
            check_row(input, first_year + year, expected[cell][year], headed);
        check_bulk(input, cell, yearly);
    }
    CHECK_FALSE(input.LoadNext());
}

TEST_CASE("TimeDataD memory loading retains the requested subset of cells", "[timedata]") {
    InputFile file(std::string("Lon Lat Year A B C\n") + local_yearly_rows);
    ScopedMapTextFile mapping(false);
    ListArray_id<Coord> gridlist;
    add_coordinate(gridlist, 2);
    add_coordinate(gridlist, 0);
    TimeDataD input;
    REQUIRE(input.Open(file.path.c_str(), gridlist));
    for (std::size_t cell : {2, 0, 2}) {
        REQUIRE(input.Load(coordinate(cell)));
        for (int year = 0; year < 3; ++year)
            check_row(input, first_year + year, expected[cell][year]);
    }
}

TEST_CASE("TimeDataD finds positive values across selected cells and years", "[timedata]") {
    InputFile file(
        "Lon Lat Year Later Zero Negative Elsewhere\n"
        "12.25 -33.75 2000 0 0 -1 0\n"
        "12.25 -33.75 2001 0 0 -2 0\n"
        "12.25 -33.75 2002 0 0 -3 0\n"
        "-120.25 45.25 2000 0 0 -4 0\n"
        "-120.25 45.25 2001 0 0 -5 0\n"
        "-120.25 45.25 2002 7 0 -6 0\n"
        "150.25 -20.25 2000 0 0 -7 8\n"
        "150.25 -20.25 2001 0 0 -8 9\n"
        "150.25 -20.25 2002 0 0 -9 10\n");
    ListArray_id<Coord> gridlist;
    add_coordinate(gridlist, 0);
    add_coordinate(gridlist, 1);
    Access access = Access::Direct;
    SECTION("direct") {}
    SECTION("mapped") { access = Access::Mapped; }
    SECTION("memory") { access = Access::Memory; }
    ScopedMapTextFile mapping(access == Access::Mapped);
    TimeDataD input;
    if (access == Access::Direct)
        REQUIRE(input.Open(file.path.c_str()));
    else
        REQUIRE(input.Open(file.path.c_str(), gridlist));
    input.CheckIfPresent(gridlist);
    char later[] = "Later", zero[] = "Zero", negative[] = "Negative", elsewhere[] = "Elsewhere";
    CHECK(input.item_has_data(later));
    CHECK_FALSE(input.item_has_data(zero));
    CHECK_FALSE(input.item_has_data(negative));
    CHECK_FALSE(input.item_has_data(elsewhere));
    REQUIRE(input.Load(coordinate(1)));
    CHECK(input.Get(2002, "Later") == Approx(7));
}

TEST_CASE("TimeDataD finds positive values in global static and yearly data", "[timedata]") {
    bool yearly = true;
    SECTION("yearly positive value occurs after the first year") {}
    SECTION("static") { yearly = false; }
    InputFile file(yearly
        ? "Year Positive Zero Negative\n2000 0 0 -1\n2001 7 0 -2\n"
        : "static Positive Zero Negative\n0 7 0 -2\n");
    TimeDataD input;
    REQUIRE(input.Open(file.path.c_str()));
    ListArray_id<Coord> gridlist;
    add_coordinate(gridlist, 0);
    input.CheckIfPresent(gridlist);
    char positive[] = "Positive", zero[] = "Zero", negative[] = "Negative";
    CHECK(input.item_has_data(positive));
    CHECK_FALSE(input.item_has_data(zero));
    CHECK_FALSE(input.item_has_data(negative));
}

TEST_CASE("TimeDataD identifies ordinary normalized fractions without altering values", "[timedata]") {
    bool global = false, yearly = true;
    SECTION("global static") { global = true; yearly = false; }
    SECTION("global yearly") { global = true; }
    SECTION("local static") { yearly = false; }
    SECTION("local yearly") {}
    for (bool normalized : {true, false}) {
        INFO("normalized=" << normalized);
        const std::string values = normalized ? "0.20 0.30 0.50\n" : "2 3 5\n";
        std::string contents;
        if (global)
            contents = yearly ? "Year A B C\n2000 " + values + "2001 " + values
                              : "static A B C\n0 " + values;
        else
            contents = yearly ? "Lon Lat Year A B C\n12.25 -33.75 2000 " + values +
                                "12.25 -33.75 2001 " + values
                              : "Lon Lat A B C\n12.25 -33.75 " + values;
        InputFile file(contents);
        TimeDataD input;
        REQUIRE(input.Open(file.path.c_str()));
        CHECK(input.NormalisedData() == normalized);
        REQUIRE(input.Load(coordinate(0)));
        check_row(input, 2000, normalized ? Row{{0.2, 0.3, 0.5}} : Row{{2, 3, 5}});
        if (yearly)
            check_row(input, 2001, normalized ? Row{{0.2, 0.3, 0.5}} : Row{{2, 3, 5}});
    }
}

TEST_CASE("TimeDataD reports decimal precision of values rather than coordinates", "[timedata]") {
    std::string contents;
    int precision = 4;
    SECTION("global static") { contents = "static A B\n0 0.10 2.3456\n"; }
    SECTION("global yearly") { contents = "Year A B\n2000 0.10 2.3456\n2001 0.20 3.000\n"; }
    SECTION("local static") { contents = "Lon Lat A B\n12.250000 -33.750000 0.10 2.3456\n"; }
    SECTION("local yearly") { contents = "Lon Lat Year A B\n12.250000 -33.750000 2000 0.10 2.3456\n"; }
    SECTION("integer values") { contents = "Year A B\n2000 1 2\n2001 3 4\n"; precision = 0; }
    InputFile file(contents);
    TimeDataD input;
    REQUIRE(input.Open(file.path.c_str()));
    CHECK(input.GetPrecision() == precision);
}

// These assert intended behavior, not today's incorrect results. They are
// hidden from the default suite and deliberately fail when selected explicitly.
// Promote them to the normal suite as the corresponding defects are resolved.
TEST_CASE("TimeDataD reads a headerless global static record", "[timedata][.timedata-known-defect]") {
    InputFile file("0 0 -2.5 1e-3\n");
    TimeDataD input(TextInput::GLOBAL_STATIC);
    REQUIRE(input.Open(file.path.c_str()));
    CHECK(input.GetNCells() == 1);
    CHECK(input.GetnYears() == 1);
    check_row(input, 2000, expected[0][0], false);
}

TEST_CASE("TimeDataD counts and reads every headerless local static row", "[timedata][.timedata-known-defect]") {
    InputFile file(local_static_rows);
    TimeDataD input(TextInput::LOCAL_STATIC);
    REQUIRE(input.Open(file.path.c_str()));
    CHECK(input.GetNCells() == 3);
    CHECK(input.GetnColumns() == 3);
    CHECK(input.GetnYears() == 1);
    for (std::size_t cell : {2, 0, 1}) {
        REQUIRE(input.Load(coordinate(cell)));
        check_row(input, 2000, expected[cell][0], false);
    }
}

TEST_CASE("TimeDataD bulk access follows memory cell selection", "[timedata][.timedata-known-defect]") {
    bool yearly = true;
    SECTION("yearly") {}
    SECTION("static") { yearly = false; }
    InputFile file(yearly ? std::string("Lon Lat Year A B C\n") + local_yearly_rows
                          : std::string("Lon Lat A B C\n") + local_static_rows);
    TimeDataD input;
    REQUIRE(open_local(input, file, Access::Memory));
    REQUIRE(input.Load(coordinate(0)));
    check_bulk(input, 0, yearly);
}

TEST_CASE("TimeDataD reports the selected coordinate for memory data", "[timedata][.timedata-known-defect]") {
    InputFile file(std::string("Lon Lat Year A B C\n") + local_yearly_rows);
    TimeDataD input;
    REQUIRE(open_local(input, file, Access::Memory));
    REQUIRE(input.Load(coordinate(0)));
    CHECK(input.GetCoord().lon == Approx(coordinates[0][0]));
    CHECK(input.GetCoord().lat == Approx(coordinates[0][1]));
}

TEST_CASE("TimeDataD reports the selected coordinate for headerless yearly data", "[timedata][.timedata-known-defect]") {
    InputFile file(local_yearly_headerless);
    TimeDataD input(TextInput::LOCAL_YEARLY);
    REQUIRE(input.Open(file.path.c_str()));
    REQUIRE(input.Load(coordinate(0)));
    CHECK(input.GetCoord().lon == Approx(coordinates[0][0]));
    CHECK(input.GetCoord().lat == Approx(coordinates[0][1]));
}
