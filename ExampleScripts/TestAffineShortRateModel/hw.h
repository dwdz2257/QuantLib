#include <ql/math/distributions/normaldistribution.hpp>
#include <ql/math/randomnumbers/boxmullergaussianrng.hpp>
#include <ql/math/randomnumbers/mt19937uniformrng.hpp>
#include <ql/math/randomnumbers/randomsequencegenerator.hpp>
#include <ql/math/randomnumbers/rngtraits.hpp>
#include <ql/methods/montecarlo/mctraits.hpp>
#include <ql/methods/montecarlo/pathgenerator.hpp>
#include <ql/models/shortrate/onefactormodels/hullwhite.hpp>
#include <ql/processes/hullwhiteprocess.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/termstructures/yield/discountcurve.hpp>
#include <cmath>
#include <chrono>
using namespace QuantLib;


void testHullWhiteModel() {
    // creat a flat forward curve

    auto act365f = Actual365Fixed();
    auto valuation_date = Date(11, Dec, 2023);
    auto reference_date = valuation_date;
    Settings::instance().evaluationDate() = valuation_date;
    double const_r = 0.03;
    auto ffc = FlatForward(reference_date, const_r, act365f);
    double a = 0.1;
    double vol = 0.015;
    auto ffc_handle = Handle<YieldTermStructure>(ext::make_shared<FlatForward>(ffc));
    auto hwmodel = HullWhite(ffc_handle, a, vol);

    auto oney_date = Date(7, Dec, 2024);
    auto oney_time = act365f.yearFraction(valuation_date, oney_date);

    std::cout << hwmodel.discount(oney_time) << "  " << exp(-const_r * oney_time) << std::endl;


    std::cout << "initial r: " << hwmodel.r0() << std::endl;

    // price bond option
    auto option_type = Option::Call;
    auto strike = 0.98;
    double tau = 1.;
    double bond_mat_short = 0.4;
    double bond_mat = 1.5;


    std::cout << "discount bond price at time zero: " << hwmodel.discountBond(0., bond_mat, const_r)
              << std::endl;
    std::cout << "discount bond price at time bond mat short: "
              << hwmodel.discountBond(bond_mat_short, bond_mat, const_r) << std::endl;
    std::cout << "forward bond price: "
              << hwmodel.discountBond(0., bond_mat, const_r) /
                     hwmodel.discountBond(0, bond_mat_short, const_r)
              << ", discount bond price mat bond short mat at time zero: "
              << hwmodel.discountBond(0, bond_mat_short, const_r) << std::endl;
    std::cout << "discount bond option price"
              << hwmodel.discountBondOption(option_type, strike, tau, bond_mat) << std::endl;

    // price bond option via mc
    auto hw_process = ext::make_shared<HullWhiteProcess>(ffc_handle, a, vol);
    size_t num_time_steps = 1750;
    size_t num_paths = 10000;
    auto gaussian_rsg = PseudoRandom::make_sequence_generator(num_time_steps, 123);
    auto hw_path_generator = SingleVariate<PseudoRandom>::path_generator_type(
        hw_process, tau, num_time_steps, gaussian_rsg, false);
    auto time_grid = hw_path_generator.timeGrid();
    double zcb_pv = 0.;
    auto before_t = std::chrono::steady_clock::now();
    for (int n = 0; n < num_paths; ++n) {
        auto path = hw_path_generator.next().value;
        double inted = 0.;
        for (int k = 0; k < num_time_steps; k++) {
            auto dt = time_grid.dt(k);
            auto sr = path.at(k);
            //std::cout << sr << std::endl;
            inted += sr * dt;
        }
        zcb_pv = (zcb_pv * n + exp(-inted)) / (n + 1);
    }
    auto after_t = std::chrono::steady_clock::now();
    std::cout << "simulation cost: " << std::chrono::duration<double>(after_t - before_t).count() << std::endl;
    std::cout << "zero coupon bond pv from mc: " << zcb_pv << ", zero coupon bond pv from model"
              << hwmodel.discountBond(0, tau, const_r) << std::endl;

    auto hwf_process = ext::make_shared<HullWhiteForwardProcess>(ffc_handle, a, vol);

    std::cout << "bond call option pv with 0 strike for notebook: "
              << hwmodel.discountBondOption(Option::Call, 0., 0.664, 0.664, 1.336) << std::endl;
    std::cout << "discount bond pv: " << hwmodel.discountBond(0, 1.336, 0.03) << std::endl;

    std::cout << "HW DONE!" << std::endl;
}


void testHullWhiteModel2() {
    std::vector<double> zeros = {-0.00217, -0.00163, -0.00119, -0.00025, 0.0006,  0.00177, 0.00302,
                                 0.00437,  0.00558,  0.00676,  0.0077,   0.00849, 0.00925, 0.00991,
                                 0.01048,  0.01087,  0.01121,  0.01152,  0.01179, 0.01204};
    auto valuation_date = Date(11, Dec, 2023);
    auto dc = Actual365Fixed();
    std::vector<Date> tenor_dates = {valuation_date};
    std::vector<double> discount_bonds = {1.};
    for (auto z : zeros) {
        tenor_dates.push_back(tenor_dates.back() + Period(365, Days));
        auto dc_factor = dc.yearFraction(valuation_date, tenor_dates.back());
        discount_bonds.push_back(pow(1 + z, -dc_factor));
        std::cout << dc_factor << ",  " << discount_bonds.back() << std::endl;
    }
    auto curve = DiscountCurve(tenor_dates, discount_bonds, dc);
    auto curve_handle = Handle<YieldTermStructure>(ext::make_shared<DiscountCurve>(curve));
    auto hw_model = HullWhite(curve_handle, 0.1, 0.016);
    std::cout <<std::setprecision(14)<< hw_model.discountBondOption(Option::Call, 0.9, 5., 6.) << std::endl;

    std::cout << "HW2 DONE!" << std::endl;


}