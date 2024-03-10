///////////////////////////////////////////////////////////////////////////////////////
/// \file outputchannel.cpp
/// \brief Classes for formatting and printing output from the model
///
/// \author Joe Siltberg
/// $Date$
///
/// This Source Code Form is subject to the terms of the Mozilla Public
/// License, v. 2.0. If a copy of the MPL was not distributed with this
/// file, You can obtain one at http://mozilla.org/MPL/2.0/.
///
///////////////////////////////////////////////////////////////////////////////////////

#include "config.h"
#include "guess.h"
#include "outputchannel.h"
#include <vector>

namespace GuessOutput {

ColumnDescriptor::ColumnDescriptor(const char* title, 
                                   int width, 
                                   int precision) 
		  : t(title),
			 w(width),
			 p(precision) {
}

const std::string& ColumnDescriptor::title() const {
	 return t;
}

int ColumnDescriptor::width() const {
	 return w;
}

int ColumnDescriptor::precision() const {
	 return p;
}

ColumnDescriptors::ColumnDescriptors(const std::vector<std::string>& titles,
                                     int width, 
                                     int precision) {
	 for (size_t i = 0; i < titles.size(); i++) {
		  columns.push_back(ColumnDescriptor(titles[i].c_str(), 
		                                     width, 
		                                     precision));
	 }
}

void ColumnDescriptors::operator+=(const ColumnDescriptor& col) {
	 columns.push_back(col);
}

void ColumnDescriptors::operator+=(const ColumnDescriptors& cols) {
	 columns.insert(columns.end(), cols.columns.begin(), cols.columns.end());
}

size_t ColumnDescriptors::size() const {
	 return columns.size();
}

const ColumnDescriptor& ColumnDescriptors::operator[](size_t i) const {
	 return columns[i];
}

TableDescriptor::TableDescriptor(const char* name,
                                 const ColumnDescriptors& columns)
		  : n(name),
			 cols(columns) {
}

const std::string& TableDescriptor::name() const {
	 return n;
}

const ColumnDescriptors& TableDescriptor::columns() const {
	 return cols;
}

Table::Table()
		  : identifier(-1) {
}

Table::Table(int id) 
		  : identifier(id) {
}

int Table::id() const {
	 return identifier;
}

bool Table::invalid() const {
	 return identifier == -1;
}

Table OutputChannel::create_table(const TableDescriptor& descriptor) {
	 Table table((int) table_descriptors.size());

	 table_descriptors.push_back(descriptor);
	 values.resize(table_descriptors.size());
	 return table;
}

void OutputChannel::add_value(const Table& table, double d) {
	 // do nothing for unused tables
	 if (table.invalid()) {
		  return;
	 }

	 // make sure this row isn't already full
	 const TableDescriptor& td = table_descriptors[table.id()];
	 if (values[table.id()].size() == td.columns().size()) {
		  fail("Added too many values to a row in table %s!", 
		       td.name().c_str());
	 }

	 values[table.id()].push_back(d);
}

const TableDescriptor& 
OutputChannel::get_table_descriptor(const Table& table) const {
	 return table_descriptors[table.id()];
}

const std::vector<double> 
OutputChannel::get_current_row(const Table& table) const {
	 return values[table.id()];
}

void OutputChannel::clear_current_row(const Table& table) {
	 values[table.id()].clear();
	 std::vector<double>().swap(values[table.id()]); // clear array memory
}

FileOutputChannel::FileOutputChannel(const char* out_dir,
                                     int coords_precision)
		  : output_directory(out_dir) {

	 // calculate suitable width for the coords columns,
	 // longitudes take at most 4 characters (-180) before the decimal 
	 // point, add the decimal point, coords_precision and a little margin:
	 const int LON_MAX_LEN = 4;
	 const int MARGIN = 2;
	 int coords_width = LON_MAX_LEN+1+coords_precision+MARGIN;

	 xtring str;
	 str.printf("%%%ds", coords_width);
	 coords_title_format = (char*)str;

	 str.printf("%%%d.%df", coords_width, coords_precision);
	 coords_format = (char*)str;
}

FileOutputChannel::~FileOutputChannel() {
	 for (size_t i = 0; i < files.size(); i++) {
		 if(files[i] != NULL)
			fclose(files[i]);
	 }
}

Table FileOutputChannel::create_table(const TableDescriptor& descriptor) {
	 Table table;

	 FILE* file = NULL;
	 if (descriptor.name() != "") {
		  std::string full_path = output_directory + descriptor.name();
		  file = fopen(full_path.c_str(), "w");
		  if (file == NULL) {
				fail("Could not open %s for output\n"\
				     "Close the file if it is open in another application",
				     full_path.c_str());
		  }
		  else {
				table = OutputChannel::create_table(descriptor);
				files.push_back(file);
				printed_header.push_back(false);
		  }
	 }

	 return table;
}

void FileOutputChannel::finish_row(const Table& table, 
                                   double lon, 
                                   double lat,
                                   int year) {

	 finish_row(table, lon, lat, year, -1, false);
}

void FileOutputChannel::finish_row(const Table& table, 
                                   double lon, 
                                   double lat,
                                   int year,
                                   int day) {

	 finish_row(table, lon, lat, year, day, true);
}

void FileOutputChannel::close_table(Table& table) {

	 // do nothing for unused tables
	 if (table.invalid()) {
		  return;
	 }

	 FILE* file = files[table.id()];
	 fclose(file);
	 files[table.id()] = NULL;
}

void FileOutputChannel::finish_row(const Table& table, 
                                   double lon, 
                                   double lat,
                                   int year, 
                                   int day, 
                                   bool print_day) {
	 // do nothing for unused tables
	 if (table.invalid()) {
		  return;
	 }

	 FILE* file = files[table.id()];

	 // make sure all columns have been added
	 const std::vector<double>& row = get_current_row(table);
	 const TableDescriptor& td = get_table_descriptor(table);

	if (row.size() == 0) {
		return;
	}

	 if (row.size() < td.columns().size()) {
		  fail("Too few values in a row in table %s", td.name().c_str());
	 }

	 // print the header if this is the first output for this file
	 if (!printed_header[table.id()]) {
		  // print title for coordinates and time columns
		  fprintf(file, coords_title_format.c_str(), "Lon");
		  fprintf(file, coords_title_format.c_str(), "Lat");
		  fprintf(file, "%8s", "Year");
		  if (print_day) {
				fprintf(file, "%8s", "Day");
		  }

		  // print each column title
		  int nbr_cols = (int) get_table_descriptor(table).columns().size();
		  for (int i = 0; i < nbr_cols; i++) {
				fputs(format_header(table, i), file);
		  }
		  fprintf(file, "\n");

		  printed_header[table.id()] = true;
	 }

	 // print out coordinates and time
	 fprintf(file, coords_format.c_str(), lon);
	 fprintf(file, coords_format.c_str(), lat);
	 fprintf(file, "%8d", year);
	 if (print_day) {
		  fprintf(file, "%8d", day);
	 }

	 // print out the values
	 for (size_t i = 0; i < row.size(); i++) {
		  fprintf(file, format(table, (int)i), row[i]);
	 }
	 fprintf(file, "\n");
	 fflush(file);

	 // start on a new row
	 clear_current_row(table);
}


const char* FileOutputChannel::format(const Table& table, int column) {
	 const TableDescriptor& td = get_table_descriptor(table);
	 const ColumnDescriptor& cd = td.columns()[column];

	 // NB: not thread safe
	 static char buf[100];

	 sprintf(buf, "%%%d.%df", cd.width(), cd.precision());

	 return buf;
}

const char* FileOutputChannel::format_header(const Table& table, int column) {
	 const TableDescriptor& td = get_table_descriptor(table);
	 const ColumnDescriptor& cd = td.columns()[column];

	 // NB: not thread safe
	 static char format[100];
	 static char buf[100];

	 sprintf(format, "%%%ds", cd.width());
	 sprintf(buf, format, cd.title().c_str());

	 return buf;	 
}

OutputRows::OutputRows(OutputChannel* output_channel, 
                       double longitude, 
                       double latitude, 
                       int year)
		  : out(output_channel),
			 lon(longitude),
			 lat(latitude),
			 y(year),
			 d(-1) {
}

OutputRows::OutputRows(OutputChannel* output_channel, 
							  double longitude, 
							  double latitude, 
							  int year, 
							  int day)
		  : out(output_channel),
			 lon(longitude),
			 lat(latitude),
			 y(year),
			 d(day) {
}

OutputRows::~OutputRows() {
	 for (size_t i = 0; i < used_tables.size(); i++) {
		  if (used_tables[i]) {
				if (d == -1) {
					 out->finish_row(Table((int)i), lon, lat, y);
				}
				else {
					 out->finish_row(Table((int)i), lon, lat, y, d);
				}
		  }
	 }
}

void OutputRows::add_value(const Table& table,double d) {
	 // do nothing for unused tables
	 if (table.invalid()) {
		  return;
	 }

	 int id = table.id();

	 // remember that this table has gotten a value
	 if (id >= static_cast<int>(used_tables.size())) {
		  used_tables.resize(id+1);
	 }
	 used_tables[id] = true;
	 
	 // send the value to the output channel
	 out->add_value(table, d);
}

#ifdef COMPRESS_OUTPUT

GZFileOutputChannel::~GZFileOutputChannel() {
    for (size_t i = 0; i < gzfiles.size(); i++) {
        if ( gzfiles[i] != NULL ) {
            gzclose(gzfiles[i]);
        } else {
            if (!printed_header[i]) {
                const TableDescriptor& descriptor = get_table_descriptor(i);
                dprintf("warning! no output created for table %s\n"
                    , descriptor.name().c_str());
            }
        }
    }
}

Table GZFileOutputChannel::create_table(const TableDescriptor& descriptor) {
    Table table;

    if (descriptor.name() != "") {
        table = OutputChannel::create_table(descriptor);
        // file will be created when first line is written
        gzfiles.push_back(NULL);
        printed_header.push_back(false);
    }

    return table;
}

void GZFileOutputChannel::close_table(Table& table) {

    // do nothing for unused tables
    if (table.invalid()) {
        return;
    }

    // untested - destructor always used
    gzFile gzfile = gzfiles[table.id()];
    if ( gzfile != NULL ) {
                gzclose(gzfile);
        gzfiles[table.id()] = NULL;
    }
    else {
        if (!printed_header[table.id()]) {
            const TableDescriptor& descriptor = get_table_descriptor(table);
            dprintf("warning! no output created for table %s\n", descriptor.name().c_str());
        }
    }
}

void GZFileOutputChannel::finish_row(const Table& table,
                                    double lon,
                                    double lat,
                                    int year) {

    finish_row(table, lon, lat, year, Date::MAX_YEAR_LENGTH - 1, false);
}

void GZFileOutputChannel::finish_row(const Table& table,
                                   double lon,
                                   double lat,
                                   int year,
                                   int day) {
    finish_row(table, lon, lat, year, day, true);
}
void GZFileOutputChannel::finish_row(const Table& table,
                                   double lon,
                                   double lat,
                                   int year,
                                   int day,
                                   bool print_day) {
    // do nothing for unused tables
    if (table.invalid()) {
        return;
    }

    gzFile gzfile = gzfiles[table.id()];

    // make sure all columns have been added
    const std::vector<double>& row = get_current_row(table);
    const TableDescriptor& td = get_table_descriptor(table);

    if (row.size() == 0) {
        return;
    }

    if (row.size() < td.columns().size()) {
        fail("Too few values in a row in table %s", td.name().c_str());
    }

    // print the header and open .gz file if this is the first output for this table
    if (!printed_header[table.id()]) {
            FILE* file = NULL;
        std::string full_path = output_directory + td.name() + ".hdr";
        file = fopen(full_path.c_str(), "w");
        if (file == NULL) {
            fail("Could not open %s for output\n"\
                    "Close the file if it is open in another application",
                    full_path.c_str());
        }
        // print title for coordinates and time columns
        fprintf(file, coords_title_format.c_str(), "Lon");
        fprintf(file, coords_title_format.c_str(), "Lat");
        fprintf(file, "%6s", "Year");

        if (print_day) {
            fprintf(file, "%5s", "Day");
        }

        // print each column title
        int nbr_cols = (int) get_table_descriptor(table).columns().size();
        for (int i = 0; i < nbr_cols; i++) {
            fputs(format_header(table, i), file);
        }
        fprintf(file, "\n");

        printed_header[table.id()] = true;
        fclose(file);

        // create .gz file
        full_path = output_directory + td.name() + ".gz";
        gzfile = gzopen(full_path.c_str(), "wb");
        if (gzfile == NULL) {
                fail("Could not open %s for output\n"		\
                "Close the file if it is open in another application",
                full_path.c_str());
        }
        else {
                gzfiles[table.id()] = gzfile;
        }
    }

    // use buffer to write only once to gzfile
//     xtring line(1024);
//     xtring buffer(1024);
     // char line[1024];
     // char buffer[1024];

    // write line to file
    if ( gzfile == NULL ) {
    fail("Problem with output file for table %s\n", td.name().c_str());
    }

    // print out coordinates and time
    gzprintf(gzfile, coords_format.c_str(), lon);
    gzprintf(gzfile, coords_format.c_str(), lat);
    gzprintf(gzfile, "%6d", year);
    if (print_day) {
        gzprintf(gzfile, "%5d", day);
    }

    // print out the values
    for (size_t i = 0; i < row.size(); i++) {
        gzprintf(gzfile, format(table, (int)i), row[i]);
    }
    gzprintf(gzfile, "\n");

     //     gzputs(gzfile, line);
    //gzfile->flush(); // don't flush every line since this degrades performance

    // start on a new row
    clear_current_row(table);
}

#endif // COMPRESS_OUTPUT

}
