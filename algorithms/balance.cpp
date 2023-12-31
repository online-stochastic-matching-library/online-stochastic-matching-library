// Apply balance to obtain weight of each offline vertex
double fill_water(vector<double> &level, const double water)
{
    double l = 0, r = water, eps = 1e-8;
    while ((r - l) > eps)
    {
        double mid = (l + r) / 2, tt = 0;
        for (auto item : level)
            if (item < mid)
                tt += mid - item;
        if (tt >= water)
            r = mid;
        else
            l = mid;
    }
    return (l + r) / 2;
}

// Match by sampling without replacement
// with weight computed by balance
vector<int> graph::balance_swor()
{
    vector<double> currentLevel(onSize + offSize, 0);
    vector<bool> selected(onSize + offSize, false);
    vector<int> res(realSize, -1);

    for (int i = 0; i < realSize; i++)
    {
        vector<double> level;

        for (int j : adj[types[i]])
            level.push_back(currentLevel[j]);

        double newLevel = fill_water(level, 1);

        double mass = 0, chosen = 0;
        vector<pair<int, double>> validMass;
        for (int j : adj[types[i]])
            if (not selected[j])
                mass += max(newLevel - currentLevel[j], 0.0);

        std::uniform_real_distribution<double> dist(0, mass);
        double sample = dist(rng);

        for (int j : adj[types[i]])
            if (not selected[j])
            {
                chosen += max(newLevel - currentLevel[j], 0.0);
                if (chosen >= sample)
                {
                    selected[j] = true;
                    res[i] = j;
                    break;
                }
            }
        for (int j : adj[types[i]])
            currentLevel[j] = max(newLevel, currentLevel[j]);
    }
    return res;
}

// Match by OCS, Huang et al. (2020)
// with weight computed by balance
vector<int> graph::balance_ocs()
{
    vector<double> currentLevel(onSize + offSize, 0);
    vector<bool> selected(onSize + offSize, false);
    vector<int> res(realSize, -1);
    auto w = [](double y)
    {
        double c = (4 - 2 * sqrt(3)) / 3;
        return exp(1.0 * y + y * y / 2.0 + c * y * y * y);
    };
    for (int i = 0; i < realSize; i++)
    {
        vector<double> level;

        for (int j : adj[types[i]])
            level.push_back(currentLevel[j]);

        double newLevel = fill_water(level, 1);

        double mass = 0, chosen = 0;
        for (int j : adj[types[i]])
            if (not selected[j])
                mass += max((newLevel - currentLevel[j]), 0.0) * w(currentLevel[j]);

        std::uniform_real_distribution<double> dist(0, mass);
        double sample = dist(rng);

        for (int j : adj[types[i]])
            if (not selected[j])
            {
                chosen += max(newLevel - currentLevel[j], 0.0) * w(currentLevel[j]);
                if (chosen >= sample)
                {
                    selected[j] = true;
                    res[i] = j;
                    break;
                }
            }
        for (int j : adj[types[i]])
            currentLevel[j] = max(newLevel, currentLevel[j]);
    }
    return res;
}
