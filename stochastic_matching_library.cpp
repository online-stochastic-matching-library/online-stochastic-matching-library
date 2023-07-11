// Main program 

#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>
#include <queue>
#include <map>
#include <string>
#include <algorithm>
#include <set>
#include "glpk.h" // For Brubach et al. (2016)

using namespace std;
mt19937 rng(random_device{}());

#include "graph.h"
#include "flow_graph.h"
#include "cycle_break_graph.h"
#include "decomposite_graph.h"
#include "natural_lp.h"
#include "read_file.cpp"
#include "algorithms/algorithms.h"

// Return size of matching
// -1 represents for not matched
int match_size(const vector<int> &res)
{
    int match = 0;
    for (int i = 0; i < (int)res.size(); i++)
        if (res[i] != -1)
            match++;
    return match;
}

//Compute mean and standard deviation
pair<double, double> compute_mean_std(const vector<double> &res, bool average = true)
{
    double mean = 0;
    for (auto item : res)
        mean += item;
    mean /= res.size();
    
    double minV = res[0];
    for (auto item : res)
        minV = min(minV, item);

    double std = 0;
    for (auto item : res)
        std += (item - mean) * (item - mean);
    
    std = sqrt(std / (res.size() - 1));
    
    if (average)
        return make_pair(mean, std);
    else
        return make_pair(minV, 0);
}


//Store the results
struct resAlg
{
    // Algorithm name
    string name;
    
    // Results of runs on one type graph (sampling online vertices)
    vector<double> resRun;
    
    // Results of runs on one dataset (sampling type graphs)
    vector<pair<double, double>> resSample;
    
    // Results of runs on all datasets
    vector<pair<double, double>> resDataset;

    resAlg(string s)
    {
        name = s;
        resRun = {};
        resSample = {};
        resDataset = {};
    }
    
    // Add result of one run
    void add_run(double item)
    {
        resRun.push_back(item);
    }
    
    // Summarize runs on one type graph (sampling online vertices)
    void summary_run(double base = 1)
    {
        auto mean_std = compute_mean_std(resRun);
        resSample.push_back(make_pair(mean_std.first / base, mean_std.second / base));
        resRun.clear();
    }

    // Summarize runs on one dataset (sampling type graphs)
    void summary_sample(bool stochastic = true)
    {
        if (resSample.size() == 1)
        {
            resDataset.push_back(resSample[0]);
            resSample.clear();
            return;
        }
        
        vector<double> res;
        for (auto item : resSample)
            res.push_back(item.first);
        
        resDataset.push_back(compute_mean_std(res, stochastic));
        resSample.clear();
    }

}   OPT("OPT"),
    stochasticSWOR("StochasticSWOR"),
    regGreedy("RegularizedGreedy"),
    ranking("Ranking"),
    balanceSWOR("Balance-SWOR"),
    balanceOCS("Balance-OCS"),
    minDegree("MinDegree"),
    feldmanMMM("FeldManEtAl"),
    bahmaniKapralov("BahmaniKapralov"),
    heaupler("HeauplerEtAl"),
    manshadiGS("ManshadiEtAl"),
    jailletLu("JailletLu"),
    jailletLuNonInt("JailletLuNonInt"),
    brubachSSX("BrubachEtAl"),
    correlated("CorrelatedSampling"),
    topHalf("TopHalfSampling"),
    poissonOCS("PoissonOCS");   //All algorithms


//Presentation order in output
vector<resAlg*> resPointer = {
    &OPT,
    &regGreedy,
    &stochasticSWOR,
    &poissonOCS,
    &minDegree,
    &balanceSWOR,
    &balanceOCS,
    &ranking,
    &heaupler,
    &topHalf,
    &correlated,
    &manshadiGS,
    &brubachSSX,
    &jailletLu,
    &jailletLuNonInt,
    &bahmaniKapralov,
    &feldmanMMM
};

vector <string> datasetName;

// Save results to directory
void save_results_to_files(string directory)
{
    ofstream fileResMean, fileResStd;
    cout << "Save results into file " << directory + "/" + "(resMean.txt,resStd.txt)" << endl;
    fileResMean.open(directory + "/" + "resMean.txt");
    fileResStd.open(directory + "/" + "resStd.txt");
    
    fileResMean << fixed << setprecision(3);
    fileResStd << fixed << setprecision(3);
    
    fileResMean << "Algorithm";
    for (auto item : datasetName)
        fileResMean << " " << item;
    fileResMean << endl;
    
    fileResStd << "Algorithm";
    for (auto item : datasetName)
        fileResStd << " " << item;
    fileResStd << endl;
    
    for (auto i : resPointer)
        if (i != &OPT)
        {
            fileResMean << (*i).name;
            fileResStd << (*i).name;
            
            for (auto j : (*i). resDataset)
            {
                fileResMean << " " << j.first;
                fileResStd << " " << j.second;
            }
            
            fileResMean << endl;
            fileResStd << endl;
            
        }
    fileResMean.close();
    fileResStd.close();

    cout << "Output Results Done!" << endl;
}

// Apply numSample runs of all algorithms on a type graph (sampling online vertices)
// Algorithms may use natural LP solution as parameter, or matching probability matrix simulated by Monte-Carlo
// Extremely slow to compute natural LP, so only set useNatural to true for small graphs
void run_on_graph(graph &g, int numSample, bool useNatural = false)
{
    // Preprocessing
    int realSize = g.online_size();
    
    map<pair<int, int>, double> typeProb  = g.optimal_matching_prob(numSample, realSize);

    natural_lp lp(g.get_adj(),g.online_size());
    map<pair<int, int>, double> naturalProb;
    if (useNatural)
        naturalProb = lp.solve_lp();
    
    vector<double> offMass;
    if (useNatural)
        offMass = g.poisson_offline_mass(naturalProb);
    else
        offMass = g.poisson_offline_mass(typeProb);
    

    vector<int> blueF, redF;
    tie(blueF, redF) = g.feldman_et_al_color();

    vector<int> blueB, redB;
    tie(blueB, redB) = g.bahmani_kapralov_color();
    
    
    vector<vector<int>> jlList = g.jaillet_lu_list();
    map<pair<int, int>, double> jlProb = g.jaillet_lu_non_integral();
    
    map<pair<int, int>, double> brubachLp = g.brubach_et_al_lp();
    vector<vector<pair<int, double>>> brubachSSXH = g.brubach_et_al_h(brubachLp);
    
    vector<int> heauplerM1, heauplerM2;
    vector<pair<int, int>> heauplerM3;
    tie(heauplerM1, heauplerM2, heauplerM3) = g.haeupler_et_al_advice(brubachLp);

    for (int i = 0; i < numSample; i++)
    {
        g.realize(realSize);
        
        OPT.add_run(match_size(g.maximum_matching()));
        
        if (useNatural)
        {
            stochasticSWOR.add_run(match_size(g.sampling_without_replacement(naturalProb)));
            regGreedy.add_run(match_size(g.regularized_greedy(naturalProb)));
            poissonOCS.add_run(match_size(g.poisson_ocs(offMass, naturalProb)));
            topHalf.add_run(match_size(g.top_half_sampling(naturalProb)));
            correlated.add_run(match_size(g.correlated_sampling(naturalProb)));
        }
        else
        {
            stochasticSWOR.add_run(match_size(g.sampling_without_replacement(typeProb)));
            regGreedy.add_run(match_size(g.regularized_greedy(typeProb)));
            poissonOCS.add_run(match_size(g.poisson_ocs(offMass, typeProb)));
            topHalf.add_run(match_size(g.top_half_sampling(typeProb)));
            correlated.add_run(match_size(g.correlated_sampling(typeProb)));
        }
        
        ranking.add_run(match_size(g.ranking()));
        balanceSWOR.add_run(match_size(g.balance_swor()));
        balanceOCS.add_run(match_size(g.balance_ocs()));
        minDegree.add_run(match_size(g.min_degree()));
        
        feldmanMMM.add_run(match_size(g.feldman_et_al(blueF, redF)));
        bahmaniKapralov.add_run(match_size(g.bahmani_kapralov(blueB, redB)));
        heaupler.add_run(match_size(g.haeupler_et_al(heauplerM1, heauplerM2, heauplerM3)));
        manshadiGS.add_run(match_size(g.manshadi_et_al(typeProb)));
        jailletLu.add_run(match_size(g.jaillet_lu(jlList)));
        jailletLuNonInt.add_run(match_size(g.manshadi_et_al(jlProb)));
        brubachSSX.add_run(match_size(g.brubach_et_al(brubachSSXH)));
    }
    
    // Summarize runs on one type graph (sampling online vertices)
    double opt = compute_mean_std(OPT.resRun).first;
    for (auto i : resPointer)
        (*i).summary_run(opt);
}

// Apply numSample runs of algorithms on a graph in online matching
// Only four algorithms available in online matching: MinDegree, RANKING, Balance-OCS, Balance-SWOR
void run_on_non_stochastic_graph(graph &g, int numSample)
{
    // Preprocessing
    int realSize = g.online_size();
    g.realize(realSize, false);
    OPT.add_run(match_size(g.maximum_matching()));

    for (int i = 0; i < numSample; i++)
    {
        ranking.add_run(match_size(g.ranking()));
        balanceSWOR.add_run(match_size(g.balance_swor()));
        balanceOCS.add_run(match_size(g.balance_ocs()));
        minDegree.add_run(match_size(g.min_degree()));
    }
    
    // Summarize runs on one type graph (sampling online vertices)
    double opt = compute_mean_std(OPT.resRun).first;
    for (auto i : resPointer)
        (*i).summary_run(opt);
}


// Run experiments on graphs generated from file
// Extremely slow to compute natural LP, so only set useNatural to true for small graphs
void work_from_file(string name, bool useNatural = false)
{
        
    cerr << "Working on file " << name << endl;

    int numGraph = 1;
    int numSample = 10000;
    
    cerr << "Rep";
    for (int i = 1; i <= numGraph; i++)
    {
        cout << " " << i;
        graph g = generate_from_file(name, true, 0);
        run_on_graph(g, numSample, useNatural);
    }
    cerr << endl;

    for (auto j : resPointer) (*j).summary_sample();
}

// Run non-stochastic experiments on graphs generated from file
void work_from_file_non_stochastic(string name, bool useNatural = false, bool stochastic = true)
{
        
    cerr << "Working on file " << name << endl;

    int numGraph = 1000;
    int numSample = 100;
    resPointer = {&OPT, &balanceSWOR, &balanceOCS, &minDegree, &ranking};

    cerr << "Rep";
    for (int i = 1; i <= numGraph; i++)
    {
        cout << " " << i;
        graph g = generate_from_file(name, true, 0);
        run_on_non_stochastic_graph(g, numSample);
    }
    cerr << endl;

    for (auto j : resPointer) (*j).summary_sample(false);
}

int main()
{
    // Work on real-world datasets in online stochastic matching
    vector<pair<string, string>> file_name = 
    {
        make_pair("real_world/socfb-Caltech36/socfb-Caltech36.txt", "Caltech36"),
        make_pair("real_world/socfb-Reed98/socfb-Reed98.txt", "Reed98"),
        make_pair("real_world/bio-CE-GN/bio-CE-GN.txt", "CE-GN"),
        make_pair("real_world/bio-CE-PG/bio-CE-PG.txt", "CE-PG"),
        make_pair("real_world/econ-beause/econ-beause.txt", "beause"),
        make_pair("real_world/econ-mbeaflw/econ-mbeaflw.txt", "mbeaflw")
    };
    for (auto item : file_name)
    {
        work_from_file(item.first);
        datasetName.push_back(item.second);
    }
    save_results_to_files("real_world_result");
   
    
    
    // Work on real-world small datasets in online stochastic matching
    /*vector<pair<string, string>> file_name = 
    {
        make_pair("real_world_small/soc-firm-hi-tech/soc-firm-hi-tech.txt", "hi-tech"),
        make_pair("real_world_small/soc-physicians/soc-physicians.edges", "physicians"),
        make_pair("real_world_small/gent113/gent113.mtx", "gent113"),
        make_pair("real_world_small/lp_blend/lp_blend.mtx", "lp_blend")
    };
    for (auto item : file_name)
    {
        work_from_file(item.first);
        datasetName.push_back(item.second);
    }
    save_results_to_files("real_world_small_result");*/
    
        
    
    // Work on real-world datasets in online matching
    /*vector<pair<string, string>> file_name = 
    {
        make_pair("real_world/socfb-Caltech36/socfb-Caltech36.txt", "Caltech36"),
        make_pair("real_world/socfb-Reed98/socfb-Reed98.txt", "Reed98"),
        make_pair("real_world/bio-CE-GN/bio-CE-GN.txt", "CE-GN"),
        make_pair("real_world/bio-CE-PG/bio-CE-PG.txt", "CE-PG"),
        make_pair("real_world/econ-beause/econ-beause.txt", "beause"),
        make_pair("real_world/econ-mbeaflw/econ-mbeaflw.txt", "mbeaflw")
    };
    for (auto item : file_name)
    {
        work_from_file_non_stochastic(item.first, false);
        datasetName.push_back(item.second);
    }
    save_results_to_files("real_world_non_stochastic_result");*/
    return 0;
}
