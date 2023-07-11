// Match by Regularized Greedy
vector<int> graph::regularized_greedy(map<pair<int, int>, double> &typeProb)
{
    double theta = 0.4253;
    auto alpha = [&theta](double t)
    {
        return (1 - (1.0 / theta * exp(-(1 - log(1 - theta)) * (1 - t)) - (1 - log(1 - theta)) * exp(-1.0 / theta * (1 - t))) / (1.0 / theta - 1 + log(1 - theta)));
    };
    auto beta = [&theta](double t)
    {
        return (exp(-(1 - log(1 - theta)) * (1 - t)) - exp(-1.0 / theta * (1 - t))) / (1.0 / theta - 1 + log(1 - theta));
    };
    auto p = [&theta](double x)
    {
        return min(x / theta, 1.0);
    };
    vector<double> offlineMass(onSize + offSize, 0);
    vector<double> onlineMass(onSize + offSize, 0);
    for (int i = 0; i < onSize; i++)
    {
        for (int j : adj[i])
        {
            double mass = typeProb[make_pair(i, j)];
            offlineMass[j] += mass;
            onlineMass[i] += mass;
        }
    }
    vector<int> res(realSize, -1);
    vector<bool> matched(onSize + offSize, false);
    
    for (int i = 0; i < realSize; i++)
    {
        double minVal = 1e10, val;
        int index = -1;
        for (int j : adj[types[i]])
        {
            if (not matched[j])
            {
                double t = 1.0 * i / realSize;
                val = alpha(t) * offlineMass[j];
                for (int onlineTypeNeighbor : adj[j])
                {
                    val += beta(t) * (p(onlineMass[onlineTypeNeighbor]) - p(onlineMass[onlineTypeNeighbor] - typeProb[make_pair(onlineTypeNeighbor, j)]));
                }
                if (minVal > val)
                {
                    minVal = val;
                    index = j;
                }
            }
        }
        if (index != -1)
        {
            res[i] = index;
            matched[index] = true;
            offlineMass[index] = 0;
            for (int onlineTypeNeighbor: adj[index])
                onlineMass[onlineTypeNeighbor] -= typeProb[make_pair(onlineTypeNeighbor, index)];
        }
    }
    return res;
}
