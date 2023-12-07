#include <ql/termstructures/yield/discountcurve.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/termstructures/yield/zerocurve.hpp>
#include <cmath>
#include <iostream>
#include <vector>

using namespace QuantLib;

int main() {

    auto act365f = Actual365Fixed();
    auto valuation_date = Date(7, Dec, 2023);
    Settings::instance().evaluationDate() = valuation_date;
    auto reference_date = valuation_date;
    auto halfy_date = Date(7, Jun, 2024);
    auto halfy_time = act365f.yearFraction(valuation_date, halfy_date);
    auto oney_date = Date(7, Dec, 2024);
    auto oney_time = act365f.yearFraction(valuation_date, oney_date);
    auto oneyhalf_date = Date(7, Jun, 2025);
    auto oneyhalf_time = act365f.yearFraction(valuation_date, oneyhalf_date);
    std::cout << "halfy time: " << halfy_time << ", oney time: " << oney_time << std::endl;
    auto twoy_date = Date(7, Dec, 2025);
    auto twoy_time = act365f.yearFraction(valuation_date, twoy_date);
    auto threey_date = Date(7, Dec, 2026);
    auto threey_time = act365f.yearFraction(valuation_date, threey_date);

    // create a flat forward curve
    double const_r = 0.05;
    auto flat_foward_curve = FlatForward(reference_date, const_r, act365f, Continuous, Annual);
    auto oney_df = flat_foward_curve.discount(oney_date);
    auto oney_df2 = exp(-const_r * oney_time);
    std::cout << "curve discount: " << oney_df << ", manual discount" << oney_df2 << std::endl;

    auto oney_zr = flat_foward_curve.zeroRate(oney_date, act365f, Continuous, Annual, false);
    auto oney_zr2 = log(1 / oney_df2) / oney_time;
    std::cout << "oney zero rate from curve: " << oney_zr << ", oney zero rate manual: " << oney_zr2
              << std::endl;

    // create a fake discounting curve
    std::vector<Date> tenor_dates = {valuation_date, oney_date, twoy_date, threey_date};
    std::vector<double> dfs = {1., 0.98, 0.95, 0.92};
    auto discount_curve = DiscountCurve(tenor_dates, dfs, act365f);
    auto dc_oney_df = discount_curve.discount(oneyhalf_date);
    auto dc_oney_df2 = exp(log(dfs[1]) + (log(dfs[2]) - log(dfs[1])) / (twoy_time - oney_time) *
                                             (oneyhalf_time - oney_time));
    std::cout << "dc curve discount: " << dc_oney_df << ", manual discount" << dc_oney_df2
              << std::endl;


    std::cout << "valuation date zero rate: " << discount_curve.zeroRate(0, Continuous)
              << "oney time zero rate" << discount_curve.zeroRate(oney_time, Continuous)
              << std::endl;

    std::cout << "oney half zero rate: " << discount_curve.zeroRate(oneyhalf_time, Continuous)
              << "twoy time zero rate" << discount_curve.zeroRate(twoy_time, Continuous)
              << std::endl;

    std::cout << "twoy time zero rate manual: " << log(1 / dfs[2]) / twoy_time << std::endl;

    // create zero rate curve from dfs
    std::vector<double> zrs = {log(1 / dfs[1]) / oney_time, log(1 / dfs[1]) / oney_time,
                               log(1 / dfs[2]) / twoy_time, log(1 / dfs[3]) / threey_time};
    auto zero_curve = InterpolatedZeroCurve<Linear>(tenor_dates, zrs, act365f);
    for (auto d : tenor_dates) {
        std::cout << "df from zero curve: " << zero_curve.discount(d)
                  << ", df from discounting curve: " << discount_curve.discount(d) << std::endl;
    }
    for (auto d : tenor_dates) {
        std::cout << "zr from zero curve: " << zero_curve.zeroRate(d, act365f, Continuous)
                  << ", zr from discounting curve: "
                  << discount_curve.zeroRate(d, act365f, Continuous) << std::endl;
    }
    // for dates inbetween, dc and zc should give different results
    // since the log linear dc perform linear interpolation on -zrs*t whilt linear zc perform linear
    // interpolation on zrs
    std::cout << "oney half df from zero curve: " << zero_curve.discount(oneyhalf_date)
              << ", oney half df from discounting curve: " << discount_curve.discount(oneyhalf_date)
              << std::endl;

    std::cout << "oney half zr from zero curve: "
              << zero_curve.zeroRate(oneyhalf_date, act365f, Continuous)
              << ", oney half zr from discounting curve: "
              << discount_curve.zeroRate(oneyhalf_date, act365f, Continuous) << std::endl;
    std::cout << "DONE!" << std::endl;
    return 0;
}