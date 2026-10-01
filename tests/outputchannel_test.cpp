#include "catch.hpp"
#include "outputchannel.h"

using namespace GuessOutput;

namespace {

class CapturingOutputChannel : public OutputChannel {
public:
	std::vector<double> numbers;
	std::vector<std::string> text;

	void finish_row(const Table& table, double, double, int) override {
		numbers = get_current_row(table);
		text = get_current_row_text(table);
		clear_current_row(table);
	}

	void finish_row(const Table& table, double lon, double lat, int year, int) override {
		finish_row(table, lon, lat, year);
	}

	void close_table(Table&) override {}
};

}

TEST_CASE("OutputChannel owns text in mixed rows", "[output]") {
	CapturingOutputChannel channel;
	ColumnDescriptors columns;
	columns += ColumnDescriptor("id", 8, 0);
	columns += ColumnDescriptor("name", 12, 0, TEXT);
	columns += ColumnDescriptor("fraction", 12, 3);
	Table table = channel.create_table(TableDescriptor("stand", columns));

	std::string name = "forest";
	channel.add_value(table, 7.0);
	channel.add_value(table, name);
	channel.add_value(table, 0.25);
	name = "changed";
	channel.finish_row(table, 0, 0, 2000);

	REQUIRE(channel.text[1] == "forest");
	REQUIRE(channel.numbers[0] == 7.0);
	REQUIRE(channel.numbers[2] == 0.25);

	channel.add_value(table, 8.0);
	channel.add_value(table, std::string("pasture"));
	channel.add_value(table, 0.75);
	channel.finish_row(table, 0, 0, 2001);
	REQUIRE(channel.text[1] == "pasture");
}
