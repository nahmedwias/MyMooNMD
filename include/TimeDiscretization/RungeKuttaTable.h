#ifndef RUNGEKUTTATABLE_H
#define RUNGEKUTTATABLE_H

/** ***************************************************************************
*
* @name       RungeKuttaTable
* @brief      
*
* @author     Naveed Ahmed, Ulrich Wilbrandt 
* @date       30.01.2015
*
******************************************************************************/

#include<vector>
#include <memory>

class RungeKuttaTable
{
  protected:
    /// number of stages
    unsigned int n_stages;
    /// order of the original method
    int order_p;
    /// order of the embedded scheme
    int order_q;
    /// 0: explicit, 1: diagonally implicit, 2: fully implicit
    unsigned int type; 
    
    // 
    std::vector<std::vector<double> > a;
    std::vector<double> b;
    std::vector<double> c;
    std::vector<double> bh;
public:
    RungeKuttaTable(std::string time_disc);
    
    // getters
    unsigned int get_n_stages() const
    { return n_stages; }
    
    std::vector<std::vector<double> > get_a() const
    { return a; }
    
    double get_a(unsigned int i, unsigned int j) const
    { return a.at(i).at(j); }
    
    std::vector<double> get_b() const
    { return b; }
    
    double get_b(unsigned int i) const { return b.at(i); }
    
    std::vector<double> get_c() const
    { return c; }
    
    double get_c(unsigned int i) const
    { return c.at(i); }
    
    std::vector<double> get_bh() const
    { return bh; }
    
    double get_bh(unsigned int i) const
    { return bh.at(i); }
    
    bool is_embedded() const { return order_q != 0; }
    
    int get_order_p() const
    { return order_p; }
    
    int get_order_q() const
    { return order_q; }
    
    unsigned int get_type() const
    { return type; }
};

#endif // RUNGEKUTTATABLE_H
