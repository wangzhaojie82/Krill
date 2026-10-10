#include "header/DatasetLoader.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

bool appendPair(PairData& data,
                ExactStats& exact_stats,
                uint32_t flow_id,
                uint32_t element_id,
                uint64_t& total_count) {
    if (flow_id == 0 || element_id == 0) {
        return false;
    }
    data.emplace_back(flow_id, element_id);
    exact_stats.update(flow_id, element_id);
    ++total_count;
    return true;
}

std::vector<std::string> splitCsvLine(const std::string& line) {
    std::vector<std::string> fields;
    std::istringstream stream(line);
    std::string token;
    while (std::getline(stream, token, ',')) {
        fields.push_back(token);
    }
    return fields;
}

void printLoadedSummary(const std::string& name, const LoadedDataset& loaded) {
    std::cout << name << " loaded. data=" << loaded.data.size()
              << ", unique_pairs=" << loaded.exact_stats.uniquePairCount()
              << "\n" << std::endl;
}

} // namespace

LoadedDataset loadDatasetECommerceEvents() {
    const std::vector<std::string> file_paths = {
            "G:/Dataset/eCommerce Events History/Cosmetics/raw/2020-Feb.csv",
            "G:/Dataset/eCommerce Events History/Cosmetics/raw/2020-Jan.csv",
            "G:/Dataset/eCommerce Events History/Cosmetics/raw/2019-Dec.csv",
            "G:/Dataset/eCommerce Events History/Cosmetics/raw/2019-Nov.csv",
            "G:/Dataset/eCommerce Events History/Cosmetics/raw/2019-Oct.csv",
    };

    LoadedDataset loaded;
    uint64_t total_data_num = 0;

    for (const auto& file_path : file_paths) {
        std::ifstream file(file_path);
        if (!file.is_open()) {
            std::cerr << "Failed to open file: " << file_path << std::endl;
            continue;
        }

        uint64_t file_data_num = 0;
        std::string line;
        if (!std::getline(file, line)) {
            std::cerr << "File " << file_path << " is empty." << std::endl;
            continue;
        }

        while (std::getline(file, line)) {
            const auto fields = splitCsvLine(line);
            if (fields.size() <= 7) {
                continue;
            }

            try {
                const uint32_t flow_id = static_cast<uint32_t>(std::stoull(fields[2]));
                const uint32_t element_id = static_cast<uint32_t>(std::stoull(fields[7]));
                if (appendPair(loaded.data, loaded.exact_stats, flow_id, element_id, total_data_num)) {
                    ++file_data_num;
                }
            } catch (const std::exception&) {
                continue;
            }
        }
        std::cout << file_path << " loaded. data=" << file_data_num << std::endl;
    }

    printLoadedSummary("ECommerceEvents", loaded);
    return loaded;
}

LoadedDataset loadDatasetTenrec() {
    const std::vector<std::string> file_paths = {
            "G:/Dataset/RecSysDataset/Tenrec/QK-article.csv",
    };

    LoadedDataset loaded;
    uint64_t total_data_num = 0;

    for (const auto& file_path : file_paths) {
        std::ifstream file(file_path);
        if (!file.is_open()) {
            std::cerr << "Failed to open file: " << file_path << std::endl;
            continue;
        }

        uint64_t file_data_num = 0;
        std::string line;
        if (!std::getline(file, line)) {
            std::cerr << "File " << file_path << " is empty." << std::endl;
            continue;
        }

        while (std::getline(file, line)) {
            const auto fields = splitCsvLine(line);
            if (fields.size() <= 1) {
                continue;
            }

            try {
                const uint32_t flow_id = static_cast<uint32_t>(std::stoull(fields[1]));
                const uint32_t element_id = static_cast<uint32_t>(std::stoull(fields[0]));
                if (appendPair(loaded.data, loaded.exact_stats, flow_id, element_id, total_data_num)) {
                    ++file_data_num;
                }
            } catch (const std::exception&) {
                continue;
            }
        }
        std::cout << file_path << " loaded. data=" << file_data_num << std::endl;
    }

    printLoadedSummary("Tenrec", loaded);
    return loaded;
}

LoadedDataset loadDataSetStackOverflow() {
    const std::vector<std::string> file_paths = {
            "G:/Dataset/StackOverflow/sx-stackoverflow/sx-stackoverflow-c2a.txt",
    };

    LoadedDataset loaded;
    uint64_t total_data_num = 0;

    for (const auto& file_path : file_paths) {
        std::ifstream file(file_path);
        if (!file.is_open()) {
            std::cerr << "Failed to open file: " << file_path << std::endl;
            continue;
        }

        uint64_t file_data_num = 0;
        std::string line;
        while (std::getline(file, line)) {
            if (line.empty()) {
                continue;
            }

            for (char& ch : line) {
                if (ch == ',' || ch == ';' || ch == '\t') {
                    ch = ' ';
                }
            }

            std::istringstream iss(line);
            uint64_t src_id = 0;
            uint64_t dst_id = 0;
            uint64_t unix_ts = 0;
            if (!(iss >> src_id >> dst_id >> unix_ts)) {
                continue;
            }
            if (src_id == 0 || dst_id == 0) {
                continue;
            }

            const uint32_t flow_id = static_cast<uint32_t>(dst_id);
            const uint32_t element_id = static_cast<uint32_t>(src_id);
            if (appendPair(loaded.data, loaded.exact_stats, flow_id, element_id, total_data_num)) {
                ++file_data_num;
            }
        }
        std::cout << file_path << " loaded. data=" << file_data_num << std::endl;
    }

    printLoadedSummary("StackOverflow", loaded);
    return loaded;
}

