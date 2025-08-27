#pragma once

struct PairHash
{
    size_t operator()(const std::pair<int, int>& p) const
    {
        size_t h1 = std::hash<int>{}(p.first);
        size_t h2 = std::hash<int>{}(p.second);

        return h1 ^ (h2 << 1);
    }
};
