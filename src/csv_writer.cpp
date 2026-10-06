#include "monograph/csv_writer.hpp"

#include <iomanip>
#include <stdexcept>

namespace monograph {

    CsvWriter::CsvWriter(const std::string& path, const std::vector<std::string>& columns)
        : out_(path), ncols_(columns.size()) {
        if (!out_.is_open()) return;
        out_ << std::setprecision(12);
        for (std::size_t i = 0; i < columns.size(); ++i) {
            if (i) out_ << ',';
            out_ << columns[i];
        }
        out_ << '\n';
    }

    void CsvWriter::write_row(const std::vector<double>& values) {
        if (values.size() != ncols_) throw std::invalid_argument("CsvWriter: wrong number of values");
        for (std::size_t i = 0; i < values.size(); ++i) {
            if (i) out_ << ',';
            out_ << values[i];
        }
        out_ << '\n';
        ++rows_;
    }

}  // namespace monograph
