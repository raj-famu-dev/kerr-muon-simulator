#pragma once
#include <fstream>
#include <string>
#include <vector>

namespace monograph {

    // Minimal CSV writer: a header row, then rows of numbers.
    class CsvWriter {
    public:
        CsvWriter(const std::string& path, const std::vector<std::string>& columns);

        bool is_open() const { return out_.is_open(); }
        long rows_written() const { return rows_; }

        // Writes one row. The number of values must match the number of columns.
        void write_row(const std::vector<double>& values);

    private:
        std::ofstream out_;
        std::size_t ncols_;
        long rows_ = 0;
    };

}  // namespace monograph
