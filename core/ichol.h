#pragma once

template <typename INT>
struct static_graph_t
{
    std::vector<INT> u_, v_;
    std::vector<INT> in_degree_;
    std::vector<INT> first_, next_;
    std::vector<bool> valid_;
    INT curr_edge_;

    static_graph_t(const INT vert_num, const INT edge_num)
        : curr_edge_(0)
    {
        in_degree_.resize(vert_num);
        std::fill(in_degree_.begin(), in_degree_.end(), 0);

        u_.resize(edge_num);
        v_.resize(edge_num);
        valid_.resize(edge_num);
        first_.resize(vert_num);
        next_.resize(edge_num);
        std::fill(first_.begin(), first_.end(), -1);
        std::fill(next_.begin(), next_.end(), -1);
        std::fill(valid_.begin(), valid_.end(), false);
    }
    INT num_nodes() const
    {
        return first_.size();
    }
    void add_edge(const INT u, const INT v)
    { // u-->v
        u_[curr_edge_] = u;
        v_[curr_edge_] = v;
        valid_[curr_edge_] = true;
        ++in_degree_[v];

        next_[curr_edge_] = first_[u];
        first_[u] = curr_edge_;

        ++curr_edge_;
    }
    INT in_degree(const INT pid)
    {
        return in_degree_[pid];
    }
    void remove_out_edges(const INT pid)
    {
        INT e = first_[pid];
        while (e != -1)
        {
            valid_[e] = false;
            --in_degree_[v_[e]];
            e = next_[e];
        }
    }
};

template <class SpMat, typename int_type>
Eigen::Matrix<int_type, -1, 1>
group_by_color(const SpMat &A, const int_type begin, const int_type end)
{
    const int_type dim = end - begin;
    Eigen::Matrix<int_type, -1, 1> order(dim);

    std::vector<int_type> colored(dim, 0);
    std::vector<bool> interact(dim, 0);

    int_type e = 1, count = 0;
    while (count < dim)
    {
        std::fill(interact.begin(), interact.end(), 0);

        for (int_type j = begin; j < end; ++j)
        {
            if (colored[j - begin] == 0 && interact[j - begin] == false)
            {
                colored[j - begin] = e;
                order[count++] = j;

                for (typename SpMat::InnerIterator it(A, j); it; ++it)
                {
                    const int_type i = it.row();
                    if (i < begin)
                        continue;
                    if (i >= end)
                        break;
                    if (i != j)
                        interact[i - begin] = true;
                }
            }
        }

        ++e;
    }

    std::cout << "color num of block: " << e - 1 << std::endl;
    std::cout << order.transpose() << std::endl;
    return order;
}