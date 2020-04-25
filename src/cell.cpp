/*! @file cell.cpp
    @brief Implementation of Cell class
*/
#include "cell.hpp"

#include <wtl/iostr.hpp>
#include <wtl/random.hpp>

#include <type_traits>

namespace tumopp {

static_assert(std::is_nothrow_copy_constructible<Cell>{}, "");
static_assert(std::is_nothrow_move_constructible<Cell>{}, "");
Cell::param_type Cell::PARAM_;

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
namespace {

class GammaFactory {
  public:
    GammaFactory(double k) noexcept: shape_(k) {}
    std::gamma_distribution<double> operator()(double mu) {
        const double theta = std::max(mu / shape_, 0.0);
        return std::gamma_distribution<double>(shape_, theta);
    }
    void param(double k) {shape_ = k;}
  private:
    double shape_;
};

template <class URBG>
inline bool bernoulli(double p, URBG& engine) {
    // consume less URBG when p is set to 0 or 1.
    return p >= 1.0 || (p > 0.0 && wtl::generate_canonical(engine) < p);
}

class bernoulli_distribution {
  public:
    bernoulli_distribution(double p) noexcept: p_(p) {}
    template <class URBG>
    bool operator()(URBG& engine) const {
        return p_ >= 1.0 || (p_ > 0.0 && wtl::generate_canonical(engine) < p_);
    }
    void param(double p) {p_ = p;}
  private:
    double p_;
};

GammaFactory GAMMA_FACTORY(Cell::param().GAMMA_SHAPE);
bernoulli_distribution BERN_SYMMETRIC(Cell::param().PROB_SYMMETRIC_DIVISION);
bernoulli_distribution BERN_MUT_BIRTH(Cell::param().RATE_BIRTH);
bernoulli_distribution BERN_MUT_DEATH(Cell::param().RATE_DEATH);
bernoulli_distribution BERN_MUT_ALPHA(Cell::param().RATE_ALPHA);
bernoulli_distribution BERN_MUT_MIGRA(Cell::param().RATE_MIGRA);
std::normal_distribution<double> GAUSS_BIRTH(Cell::param().MEAN_BIRTH, Cell::param().SD_BIRTH);
std::normal_distribution<double> GAUSS_DEATH(Cell::param().MEAN_DEATH, Cell::param().SD_DEATH);
std::normal_distribution<double> GAUSS_ALPHA(Cell::param().MEAN_ALPHA, Cell::param().SD_ALPHA);
std::normal_distribution<double> GAUSS_MIGRA(Cell::param().MEAN_MIGRA, Cell::param().SD_MIGRA);
std::poisson_distribution<int> poisson_distribution (Cell::param().RATE_PASSENGER);
std::uniform_int_distribution<int> uniform_distribution(1000001, 50000000); // define the range of mutational space, here exome
std::uniform_int_distribution<int> uniform_distribution_birth(1, 200000);  // define the range of driver (birth rate) mutational space
std::uniform_int_distribution<int> uniform_distribution_migrate(200001, 500000);   // define the range of driver (death rate) mutational space
std::uniform_int_distribution<int> uniform_distribution_death(500001, 1000000);  // define the range of driver (migration rate) mutational space
  
}// namespace
/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

void Cell::param(const param_type& p) {
    PARAM_ = p;
    GAMMA_FACTORY.param(PARAM_.GAMMA_SHAPE);
    BERN_SYMMETRIC.param(PARAM_.PROB_SYMMETRIC_DIVISION);
    BERN_MUT_BIRTH.param(PARAM_.RATE_BIRTH);
    BERN_MUT_DEATH.param(PARAM_.RATE_DEATH);
    BERN_MUT_ALPHA.param(PARAM_.RATE_ALPHA);
    BERN_MUT_MIGRA.param(PARAM_.RATE_MIGRA);
    GAUSS_BIRTH.param(decltype(GAUSS_BIRTH)::param_type(PARAM_.MEAN_BIRTH, PARAM_.SD_BIRTH));
    GAUSS_DEATH.param(decltype(GAUSS_DEATH)::param_type(PARAM_.MEAN_DEATH, PARAM_.SD_DEATH));
    GAUSS_ALPHA.param(decltype(GAUSS_ALPHA)::param_type(PARAM_.MEAN_ALPHA, PARAM_.SD_ALPHA));
    GAUSS_MIGRA.param(decltype(GAUSS_MIGRA)::param_type(PARAM_.MEAN_MIGRA, PARAM_.SD_MIGRA));
    poisson_distribution.param(decltype(poisson_distribution)::param_type(PARAM_.RATE_PASSENGER));   //ruping
}

void Cell::differentiate(urbg_t& engine) {
    if (is_differentiated()) return;
    if (BERN_SYMMETRIC(engine)) return;
    proliferation_capacity_ = static_cast<int8_t>(PARAM_.MAX_PROLIFERATION_CAPACITY);
}

std::string Cell::mutate(urbg_t& engine, urbg_t& engine3) {
    auto oss = wtl::make_oss();
    if (BERN_MUT_BIRTH(engine)) {
        event_rates_ = std::make_shared<EventRates>(*event_rates_);
        double s = GAUSS_BIRTH(engine);
        oss << id_ << "\tbeta\t" << std::to_string(uniform_distribution_birth(engine3)) << "\t" << s << "\n";
        event_rates_->birth_rate *= (s += 1.0);
    }
    if (BERN_MUT_DEATH(engine)) {
        event_rates_ = std::make_shared<EventRates>(*event_rates_);
        double s = GAUSS_DEATH(engine);
        oss << id_ << "\tdelta\t" << std::to_string(uniform_distribution_death(engine3)) << "\t" << s << "\n";
        event_rates_->death_rate *= (s += 1.0);
    }
    if (BERN_MUT_ALPHA(engine)) {
        event_rates_ = std::make_shared<EventRates>(*event_rates_);
        double s = GAUSS_ALPHA(engine);
        oss << id_ << "\talpha\t" << std::to_string(uniform_distribution_death(engine3)) << "\t" << s << "\n";
        event_rates_->death_prob *= (s += 1.0);
    }
    if (BERN_MUT_MIGRA(engine)) {
        event_rates_ = std::make_shared<EventRates>(*event_rates_);
        double s = GAUSS_MIGRA(engine);
        oss << id_ << "\trho\t" << std::to_string(uniform_distribution_migrate(engine3)) << "\t" << s << "\n";
        event_rates_->migra_rate *= (s += 1.0);
    }
    return oss.str();
}
 
std::string Cell::mutate2(urbg_t& engine2, urbg_t& engine3) {
    auto oss = wtl::make_oss();
    auto n_passengers = poisson_distribution(engine2);
    std::string passengers = "";
    
    for(int n=0; n < n_passengers; ++n)
      passengers = passengers + std::to_string(uniform_distribution(engine3)) + ",";
    
    oss << id_ << "\t" << passengers << "\n";
    return oss.str();
}
  
std::string Cell::force_mutate(urbg_t& engine) {
    event_rates_ = std::make_shared<EventRates>(*event_rates_);
    const double s_birth = GAUSS_BIRTH(engine);
    const double s_death = GAUSS_DEATH(engine);
    const double s_alpha = GAUSS_ALPHA(engine);
    const double s_migra = GAUSS_MIGRA(engine);
    event_rates_->birth_rate *= (1.0 + s_birth);
    event_rates_->death_rate *= (1.0 + s_death);
    event_rates_->death_prob *= (1.0 + s_alpha);
    event_rates_->migra_rate *= (1.0 + s_migra);
    auto oss = wtl::make_oss();
    if (s_birth != 0.0) {oss << id_ << "\tbeta\t"  << s_birth << "\n";}
    if (s_death != 0.0) {oss << id_ << "\tdelta\t" << s_death << "\n";}
    if (s_alpha != 0.0) {oss << id_ << "\talpha\t" << s_alpha << "\n";}
    if (s_migra != 0.0) {oss << id_ << "\trho\t"   << s_migra << "\n";}
    return oss.str();
}

std::string Cell::seeding(unsigned int Ns) {
    auto oss = wtl::make_oss();
    oss << id_ << "\t" << Ns << "\n";
    return oss.str();
}

double Cell::delta_time(urbg_t& engine, const double now, const double positional_value, const bool surrounded) {
    double t_birth = std::numeric_limits<double>::infinity();
    double t_death = std::numeric_limits<double>::infinity();
    double t_migra = std::numeric_limits<double>::infinity();
    if (proliferation_capacity_ != 0) {
        double mu = 1.0;
        mu /= birth_rate();
        mu /= positional_value;
        if (!surrounded) mu -= (now - time_of_birth_);
        t_birth = GAMMA_FACTORY(mu)(engine);
    }
    if (death_rate() > 0.0) {
        std::exponential_distribution<double> exponential(death_rate());
        t_death = exponential(engine);
    }
    if (migra_rate() > 0.0) {
        std::exponential_distribution<double> exponential(migra_rate());
        t_migra = exponential(engine);
    }

    if (t_birth < t_death && t_birth < t_migra) {
        next_event_ = bernoulli(death_prob(), engine)
                      ? Event::death : Event::birth;
        return t_birth;
    } else if (t_death < t_migra) {
        next_event_ = Event::death;
        return t_death;
    } else {
        next_event_ = Event::migration;
        return t_migra;
    }
}

void Cell::set_cycle_dependent_death(urbg_t& engine, const double p) {
    //TODO: reduce redundant copy for susceptible cells
    event_rates_ = std::make_shared<EventRates>(*event_rates_);
    event_rates_->death_prob = p;
    next_event_ = bernoulli(p, engine) ? Event::death : Event::birth;
}

std::string Cell::header() {
    std::ostringstream oss;
    oss << "x\ty\tz\t"
        << "id\tancestor\t"
        << "birth\tdeath\t"
        << "omega\tNt";     //ruping
    return oss.str();
}

std::ostream& Cell::write(std::ostream& ost) const {
    return ost
        << coord_[0] << "\t" << coord_[1] << "\t" << coord_[2] << "\t"
        << id_ << "\t"
        << (ancestor_ ? ancestor_->id_ : 0u) << "\t"
        << time_of_birth_ << "\t" << time_of_death_ << "\t"
        << static_cast<int>(proliferation_capacity_) << "\t" << TsizeNow_ ;
}

std::ostream& Cell::traceback(std::ostream& ost, std::unordered_set<unsigned>* done) const {
    write(ost) << "\n";
    if (ancestor_ && done->insert(ancestor_->id_).second) {
        ancestor_->traceback(ost, done);
    }
    return ost;
}

//! Stream operator for debug print
std::ostream& operator<< (std::ostream& ost, const Cell& x) {
    return x.write(ost);
}

} // namespace tumopp
