
#include "../Database/Database.h"
#include "../Util/Util.h"

using namespace std;

const string QUERY_PREFIX = "./Experiments/";
const string RESULT_PREFIX = "./ExperimentResults/Total_time/";
extern long UTIL_current_query_time;
extern vector<double> test_scores;
extern vector<vector<unsigned int>> result_ids;
extern bool k_value_feed;
extern int benchmark_k;
extern bool any_k;
extern int any_k_step;
extern vector<long> any_k_time;
extern std::vector<int> fq_num_vec;
extern std::vector<int> fr_num_vec;
extern std::vector<int> ow_num_vec;
extern bool iterator_count_enable;
extern int FQ_NUM;
extern int FR_NUM;
extern int OW_NUM;
extern RootStrategy global_root_strategy;
extern string root_var_name;

string TopKStrategyString(TopKStrategy method) {
    if (method == TopKStrategy::DP_B)
        return "DP_B";
    if (method == TopKStrategy::DP_P)
        return "DP_P";
    if (method == TopKStrategy::K_TPM)
        return "kTPM";
    if (method == TopKStrategy::DFS_OLD)
        return "DFS_OLD";
    if (method == TopKStrategy::PTAB)
        return "PTAB";
    if (method == TopKStrategy::TAKE2)
        return "Take2";
    if (method == TopKStrategy::TAKE_ALL)
        return "TAKE_ALL";
    if (method == TopKStrategy::EAGER)
        return "Eager";
    if (method == TopKStrategy::SAE)
        return "SAE";
    if (method == TopKStrategy::RankAfterMatching)
        return "Naive";
    return "Unknown";
}

vector<string> db_names{string("bsbm"), string("DBpedia3.5.1"), string("Yago")};
vector<string> query_folders{string(QUERY_PREFIX + "bsbmQueries"),
                             string(QUERY_PREFIX + "DBpedia3.5.1Queries"),
                             string(QUERY_PREFIX + "YagoQueries")};

vector<TopKStrategy> method_vec{TopKStrategy::PTAB,
                                TopKStrategy::TAKE2,
                                TopKStrategy::K_TPM,
                                TopKStrategy::EAGER,
                                TopKStrategy::DP_B};
//                                TopKStrategy::RankAfterMatching};

vector<vector<string>> query_names{{"1", "2", "3", "4", "5", "6"},
                                   {"1", "2", "3", "4", "5", "6"},
                                   {"1", "2", "3", "4", "5", "6"}};

constexpr int avg_times = 3;

void TestQuery(Database &_db, ofstream &result_file_stream, const string &db_name, const string &query_folder, const string &query_name, bool warm_up) {
    string query_file = query_folder + "/q" + query_name + ".sql";
    auto query = Util::getQueryFromFile(query_file.c_str());
    if (query.empty()) {
        cout << "query file: " << query_file << " not exists." << endl;
        return;
    }
//    result_file_stream << "q" + query_name << ".sql" << endl;

    ResultSet _rs;
    FILE *ofp = stdout;
    string msg;


    if (warm_up) { // warm up
        for (auto method: method_vec) {
            cout << TopKStrategyString(method) << endl;
            if (method == TopKStrategy::RankAfterMatching) continue;
            int ret = _db.query(query, _rs, ofp, true, false, nullptr, method);
        }
    }

    for (unsigned index = 0; index < method_vec.size(); ++index) {
      auto method = method_vec[index];
        string method_name = TopKStrategyString(method);
        cout << "method name: " << method_name << endl;
        long total_time = 0;
        for (unsigned int t = 0; t < avg_times; ++t) {
            int ret = _db.query(query, _rs, ofp, true, false, nullptr, method);
            total_time += UTIL_current_query_time;
        }
        result_file_stream << (total_time * 1.0) / (avg_times * 1.0) << (index == method_vec.size()-1 ? "\n" : ",");
    }
    cout << endl;
}

int main(int argc, char *argv[]) {
    k_value_feed = false;
    Util util;
    cout << "Test Total Time..." << endl;
    if (!util.dir_exist(RESULT_PREFIX))
        util.create_dir(RESULT_PREFIX);

    ofstream result_file_stream(RESULT_PREFIX + "total_time.txt");
    for (int i = 0; i < db_names.size(); ++i) {
        string db_name = db_names[i];
        cout << "DB name: " << db_name << endl;
        result_file_stream << db_name << endl;
        Database _db(db_name);
        _db.load();
        for (unsigned int q_i = 0; q_i < query_names[i].size(); q_i++) {
            cout << "query: " << query_folders[i] << "  " << query_names[i][q_i] << "...." << endl;
            TestQuery(_db, result_file_stream, db_names[i], query_folders[i], query_names[i][q_i], true);
        }
        result_file_stream << endl;
    }

    return 0;
}
