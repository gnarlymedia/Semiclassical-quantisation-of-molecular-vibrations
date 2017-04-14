#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <cpgplot.h>
#define CURRENT_COMP_MORSE "morse"
#define TESTING "testing"

#ifndef EXIT_FAILURE
    #define EXIT_FAILURE 1
#endif
#ifndef EXIT_SUCCESS
    #define EXIT_SUCCESS 0
#endif

// general variables
double pi = 3.1415926;
double prec;
int nsteps;
char * current_comp;
char * testing;
char * print_plugin_values = "false";

// molecular vibration variables
double e_n_global, x_min_mol_vib, a_mol_vib, x_in, x_out;
int n_mol_vib;
double v_gamma = 33.6567;
double V_0_mol_vib = 4.747;
double r_min_mol_vib = 0.74166;
int numberOfNValues = 15;
double EnergyLevels_mol_vib[15] = {
	-4.477,
	-3.962,
	-3.475,
	-3.017,
	-2.587,
	-2.185,
	-1.811,
	-1.466,
	-1.151,
	-0.867,
	-0.615,
	-0.400,
	-0.225,
	-0.094,
	-0.017
};
float valuesOfN[15] = {};
float EnergyLevels_comp[15] = {};

int isDivisible(int num, int denom, int debug)
{
        int frac = 0;
        int res = 0;
        frac = num / denom;
        if (frac * denom == num)
        {
                res = 1;
        }
	if (1 == debug)
	{
        	printf("num: %d\n", num);
        	printf("denom: %d\n", denom);
        	printf("res: %d\n", res);
	}
        return res;
}

double rect_int(double low, double high, double (* func)(double))
{
	int i = 0;
	double sum = 0.0;
	double h = (high - low)/(double) nsteps;
	double x = low;
	sum = func(x);
	for (i = 1; i < nsteps; i = i + 1)
	{
		x = x + h;
		sum = sum + func(x);
	}

	sum = sum + func(x);
	return sum * h;
}

double trap_int(double low, double high, double (* func)(double))
{
	int i;
	double sum = 0.0;
	double h = (high - low)/(double) nsteps;
	double x = low;

	for (i = 1; i <= nsteps; i = i + 1)
	{
		if (i == 1 || i == nsteps)
		{
			sum = sum + func(x);
		}
        else
        {
			sum = sum + 2.0 * func(x);
		}

		x = x + h;
	}

	return sum * h * 0.5;
}

double simp_int(double low, double high, double (* func)(double))
{
	int i;
	double sum = 0.0;
	double h = (high - low)/(double) nsteps;
	double x = low;

	for (i = 1; i <= nsteps; i = i + 1)
	{
		if (i == 1 || i == nsteps)
		{
			sum = sum + func(x);
		} else {
			if (isDivisible(i, 2, 0))
			{
				sum = sum + 4.0 * func(x);
			} else
			{
				sum = sum + 2.0 * func(x);
			}
		}

		x = x + h;
	}

	return sum * h / 3.0;
}

double ln_function(double x)
{
    if (x == 0)
    {
        return 1;
    } else
    {
        return log(1.0 + x)/x;
    }
}

double parab_function(double x)
{
    return x * x;
}

double quadratic_function(double x)
{
    return 4 * (x - 1) * (x - 2);
}

double mol_vibr_quadratic_integrand_function(double x)
{
    double quadratic_value = quadratic_function(x);

    if (e_n_global - quadratic_value > 0) {
        return sqrt(e_n_global - quadratic_value);
    }
    else
    {
        return 0.0;
    }
}

double mol_vibr_function_quadratic(double e_n)
{
    e_n_global = e_n;

    double sqrt_factor = sqrt(1 + e_n);
    x_in = (3.0 - sqrt_factor) / 2;
    x_out = (3.0 + sqrt_factor) / 2;
    double gamma_const = 1.0;

    double simp_int_eval = simp_int(x_in, x_out, mol_vibr_quadratic_integrand_function);

    return gamma_const * simp_int_eval - ((double)n_mol_vib + 0.5) * pi;
}

double morse_potential_function(double x)
{
    double exp_arg = (-1.0 * x) + (r_min_mol_vib / a_mol_vib);
    double exp_term = exp(exp_arg);
    double morse_pot_val = pow((1.0 - exp_term), 2.0) - 1;
    return morse_pot_val;
}

double mol_vibr_morse_integrand_function(double x)
{
    double morse_difference = e_n_global - morse_potential_function(x);
    double return_val;

    if (morse_difference > 0.0) {
        return_val = sqrt(morse_difference);
    }
    else
    {
        return_val = 0.0;
    }

    return return_val;
}

double mol_vibr_function_morse(double e_n)
{
    // update global value of e_n with local one
    e_n_global = e_n;

    double sqrt_argument = e_n + 1;

    if (sqrt_argument < 0) {
        sqrt_argument = 0.0;
    }

    double sqrt_factor = sqrt(sqrt_argument);
    double ln_argument_x_in = 1 - sqrt_factor;

    if (ln_argument_x_in < 0)
    {
        ln_argument_x_in = 0.0;
    }

//    printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\bE_n: %lf", V_0_mol_vib * e_n);

    x_in = x_min_mol_vib - log(ln_argument_x_in);

    double ln_argument_x_out = 1 + sqrt_factor;

    if (ln_argument_x_out < 0)
    {
        ln_argument_x_out = 0;
    }

    x_out = x_min_mol_vib - log(ln_argument_x_out);

    if (x_in > x_out) {
        // swap them if they're the wrong way around
        double x_temp = x_in;
        x_in = x_out;
        x_out = x_temp;
    }

    double gamma_const = v_gamma * a_mol_vib;
    double simp_int_eval = simp_int(x_in, x_out, mol_vibr_morse_integrand_function);
    double mol_vibr_func_morse_value = gamma_const * simp_int_eval - ((double)n_mol_vib + 0.5) * pi;

    // printf("e_n: %lf\n", e_n);
    // printf("mol_vibr_func_morse_value: %lf\n", mol_vibr_func_morse_value);

    return mol_vibr_func_morse_value;
}

double integral_parab_function_rect(double x)
{
    // the function to find the root of
	return rect_int(0.0, x, parab_function) - x;
}

double integral_parab_function_trap(double x)
{
    // the function to find the root of
	return trap_int(0.0, x, parab_function) - x;
}

double integral_parab_function_simp(double x)
{
    // the function to find the root of
	return simp_int(0.0, x, parab_function) - x;
}

double find_root_false_pos(double x_current, double x_prev, double prec, double (* func)(double))
{
    // test the starting values lie on either side of the root
    if (func(x_current) * func(x_prev) < 0)
    {
        // continue while the absolute value of the function at the
        // current x value is greater than the desired precision
        while (fabs(func(x_current)) > prec)
        {
            double f_x_prev = func(x_prev);
            double f_x_current = func(x_current);

            double x_next = x_current - ((f_x_current * (x_current - x_prev))/(f_x_current - f_x_prev));

            if (func(x_current) * func(x_next) < 0)
            {
                // they are of opposite sign - use x_current and x_next
                x_prev = x_current;
                x_current = x_next;
            }
            else
            {
                // they are of opposite sign - use x_prev and x_next
                x_current = x_next;
            }
        }
    }

    return x_current;
}

double find_root(double guess, double prec, double fstep, double (* func)(double))
{
    // first find values of x for which the function straddles the x-axis
	int i = 0;
	double xprev, xn, xnext;
	double yprev, yn, ynext;
	double cstep = fstep;
	int iLoop = 1;
	xprev = guess;
	xn = guess + cstep;
	yprev = func(xprev);
	yn = func(xn);
	for (i = 0; (i < nsteps) && (yprev * yn > 0.0); i = i + 1)
    {
		if (fabs(yn) > fabs(yprev))
		{
			cstep = -cstep/2.0;
		}

		xprev= xn;
		xn = xn + cstep;
		yn = func(xn);
		yprev = func(xprev);
	}

    // use these values to find root
	return find_root_false_pos(xn, xprev, prec, func);
}

void evaluate_ln_integral()
{
    // Evaluate ln integral in Equation 4.13
    prec = 1.0E-6;
    double low_val = 0.0;
    double high_val = 1.0;
    double diff_val = high_val - low_val;
    double actual_int = pi * pi / 12.0;
    printf("Actual: integral = %lf\n", actual_int);

    printf("\n\n\nEvaluating the ln integral: \n");
    printf("------------------------\n");
    nsteps = 100000;
    printf("Number of steps: %d\n", nsteps);
    printf("Rect: integral = %lf\n", rect_int(low_val, high_val, ln_function));

    printf("\n");
    printf("Comparing Trapezoidal and Simpson rules\n\n");

    printf("\n   Step size\t|  Trap Int\t| Abs Trap Diff\t|  %% Trap Diff\t|   Simp Int\t| Abs Simp Diff\t|   %% Simp Diff\n");
    printf("------------------------------------------------------------------------------------------------\n");

    double trap_val, trap_diff, simp_val, simp_diff, step_size;
    for (nsteps = 1; nsteps <= 1000000; nsteps = nsteps * 10){
        trap_val = trap_int(low_val, high_val, ln_function);
        trap_diff = trap_val - actual_int;
        simp_val = simp_int(low_val, high_val, ln_function);
        simp_diff = simp_val - actual_int;
        step_size = diff_val / nsteps;
        printf("   %lf\t|  %lf\t|   %lf\t|  %lf %%\t|   %lf\t|  %lf\t|  %lf %%\n", step_size, trap_val, trap_diff, trap_diff * 100 / actual_int, simp_val, simp_diff, simp_diff * 100 / actual_int);
    }
}

void root_finding_using_false_position()
{
    printf("\n------------------------\nRoot Finding Using False Position\n");
    double analytic_answer = sqrt(3);
    printf("Analytic answer for root: %lf\n\n", analytic_answer);

    printf("\n Starting Pos\t| Step Size\t|   Tolerance\t| Divisions in integral\t| Simp Root\t| Abs Diff\t| Perc Diff\n");
    printf("-----------------------------------------------------------------------------------------------------------------\n");
    int guess;
    double found_root_rect, found_root_trap, found_root_simp, rf_diff_simp, rf_perc_diff_simp, root_finding_step_size;
    double found_root_simp_best = 0.0;
    double rf_diff_simp_best = 0.0;
    double rf_perc_diff_simp_best = 100.0;
    int guess_best = 0;
    double root_finding_step_size_best = 0.0;
    double prec_best = 0.0;
    int nsteps_best = 100;
    for (guess = 2; guess <= 3; guess = guess + 1)
    {
        for (root_finding_step_size = 0.1; root_finding_step_size >= 0.01; root_finding_step_size = root_finding_step_size / 10.0)
        {
            for (prec = 0.01; prec >= 1E-3; prec = prec / 10.0)
            {
                for (nsteps = 1E5; nsteps <= 1E6; nsteps = nsteps * 10)
                {
                    found_root_simp = find_root((double)guess, prec, root_finding_step_size, integral_parab_function_simp);
                    rf_diff_simp = found_root_simp - analytic_answer;
                    rf_perc_diff_simp = rf_diff_simp * 100 / analytic_answer;

                    printf("      %d \t|   %lf\t|   %lf\t| \t%d\t\t| %lf\t| %lf\t| %lf %%\n", guess, root_finding_step_size, prec, nsteps, found_root_simp, rf_diff_simp, rf_perc_diff_simp);

                    if (fabs(rf_perc_diff_simp) < fabs(rf_perc_diff_simp_best))
                    {
                        guess_best = guess;
                        root_finding_step_size_best = root_finding_step_size;
                        prec_best = prec;
                        nsteps_best = nsteps;
                        found_root_simp_best = found_root_simp;
                        rf_diff_simp_best = rf_diff_simp;
                        rf_perc_diff_simp_best = rf_perc_diff_simp;
                    }
                }
            }
        }
    }
    printf("\n----------------------------------------------------------------------------\nSmallest percentage difference:\n");
    printf("      %d \t|   %lf\t|   %lf\t| \t%d\t\t| %lf\t| %lf\t| %lf %%\n", guess_best, root_finding_step_size_best, prec_best, nsteps_best, found_root_simp_best, rf_diff_simp_best, rf_perc_diff_simp_best);

    // rectangular
    found_root_rect = find_root((double)guess_best, prec_best, root_finding_step_size_best, integral_parab_function_rect);
    double found_root_rect_diff = found_root_rect - analytic_answer;
    printf("\n\nCompare with root found using rectangular method with same values:\n");
    printf("      %d \t|   %lf\t|   %lf\t| \t%d\t\t| %lf\t| %lf\t| %lf %%\n", guess_best, root_finding_step_size_best, prec_best, nsteps_best, found_root_rect, found_root_rect_diff, found_root_rect_diff * 100 / analytic_answer);

    // trapezoidal
    found_root_trap = find_root((double)guess_best, prec_best, root_finding_step_size_best, integral_parab_function_trap);
    double found_root_trap_diff = found_root_trap - analytic_answer;
    printf("\n\nCompare with root found using trapezoidal method with same values:\n");
    printf("      %d \t|   %lf\t|   %lf\t| \t%d\t\t| %lf\t| %lf\t| %lf %%\n", guess_best, root_finding_step_size_best, prec_best, nsteps_best, found_root_trap, found_root_trap_diff, found_root_trap_diff * 100 / analytic_answer);

    // check the roots we've found
    printf("\n-----------------------\nCheck real roots\n");
    double check_real_root_rect = integral_parab_function_rect(found_root_rect);
    printf("Check_real_root_rect:  %lf\n", check_real_root_rect);

    double check_real_root_trap = integral_parab_function_trap(found_root_trap);
    printf("Check_real_root_trap:  %lf\n", check_real_root_trap);

    double check_real_root_simp = integral_parab_function_simp(found_root_simp);
    printf("Check_real_root_simp:  %lf\n", check_real_root_simp);

    double check_real_root_analytic = integral_parab_function_simp(analytic_answer);
    printf("Check_real_root_analytic:  %lf\n", check_real_root_analytic);
}

void molecular_vibration()
{
    // mol vibration
    printf("\n------------------------\nMolecular vibration\n");

    // quadratic
    printf("\nQuadratic potential\n");
    prec = 1.0E-6;
    nsteps = 10000;
    x_in = x_out = 0.0;

    printf("\n n | Tolerance\t| Divisions in integral\t| Found Root\t| Abs Diff\t| Perc Diff\n");
    printf("------------------------------------------------------------------------------------------------------\n");
    double mol_vibr_quad_analytic_val, found_root_mol_vib_quad, mol_vibr_quad_diff_val;
    for (n_mol_vib = 0; n_mol_vib < 4; n_mol_vib++)
    {
        mol_vibr_quad_analytic_val = 4.0 * n_mol_vib + 1.0;

        for (prec = 0.01; prec >= 1E-3; prec = prec / 10.0)
        {
            for (nsteps = 1E4; nsteps <= 1E5; nsteps = nsteps * 10)
            {
                found_root_mol_vib_quad = find_root(1.0, prec, 0.1, mol_vibr_function_quadratic);
                mol_vibr_quad_diff_val = found_root_mol_vib_quad - mol_vibr_quad_analytic_val;
                printf(" %d |  %lf \t|   %d\t\t|   %lf\t| %lf\t| %lf %%\n", n_mol_vib, prec, nsteps, found_root_mol_vib_quad, mol_vibr_quad_diff_val, mol_vibr_quad_diff_val * 100 / mol_vibr_quad_analytic_val);
            }
        }
    }

    // Morse
    printf("\nMorse potential\n");
    nsteps = 10000;
//    x_in = x_out = 0.0;
    current_comp = CURRENT_COMP_MORSE;
    a_mol_vib = 0.5155;
//   a_mol_vib = 1.0;
    x_min_mol_vib = r_min_mol_vib / a_mol_vib;

    printf("\n   n\t   Comp E_n\t   Exper E_n\t  Abs Diff\t  Percent Diff\n");
    printf("------------------------------------------------------------------------\n");
    //    n_mol_vib = 0;
    double found_root_mol_vib_morse;
    for (n_mol_vib = 0; n_mol_vib < numberOfNValues; n_mol_vib++)
    {
        found_root_mol_vib_morse = V_0_mol_vib * find_root(-1.0, prec, 0.1, mol_vibr_function_morse);
        double difference = found_root_mol_vib_morse - EnergyLevels_mol_vib[n_mol_vib];
        double perc_diff = difference * 100.0 / EnergyLevels_mol_vib[n_mol_vib];
        printf("   %d\t    %lf\t   %lf\t  %lf\t  %lf %%\n", n_mol_vib, found_root_mol_vib_morse, EnergyLevels_mol_vib[n_mol_vib], difference, perc_diff);

        // set comp values in array to use for plotting
        valuesOfN[n_mol_vib] = n_mol_vib;
        EnergyLevels_comp[n_mol_vib] = found_root_mol_vib_morse;
    }

    // Plotting
    if (1 != cpgbeg(0, "?", 1, 1))
//    if (1 != cpgbeg(0, "proj_1_plot.ps/VCPS", 1, 1))
//    if (1 != cpgbeg(0, "/XWINDOW", 1, 1))
    {
        exit(EXIT_FAILURE);
    }

    double x_min = 0.0;
    double x_max = 6.0;
    double y_min = -5.0;
    double y_max = 10.0;

    // set up environment
    cpgenv(x_min, x_max, y_min, y_max, 0, 2);

    cpgbbuf();

    // set colour to black
    cpgsci(1);

    // labels
    cpglab("x", " E_n ", "Experimental vibrational energy levels of the Hydrogen molecule");

    int t,u,v;
    int plotting_resolution = 60;

    // plot morse potential
    // set colour to green
    cpgsci(3);
    float x_morse[plotting_resolution], y_morse[plotting_resolution];
    int s = sizeof(x_morse) / sizeof(x_morse[0]);
    for (t = 0; t < plotting_resolution; t++)
    {
        x_morse[t] = 0.1 * t;
        y_morse[t] = V_0_mol_vib * morse_potential_function((double)x_morse[t]);
    }
    cpgline(s, x_morse, y_morse);

    // plot lines of constant energy E_n - experimental
    // set colour to blue
    cpgsci(4);
    float x_en[plotting_resolution], y_en[plotting_resolution];
    for (u = 0; u < numberOfNValues; u++)
    {
        for (v = 0; v < plotting_resolution; v++)
        {
                x_en[v] = 0.1 * v;
                y_en[v] = EnergyLevels_mol_vib[u];
        }
        int y = sizeof(x_en) / sizeof(x_en[0]);
        cpgline(y, x_en, y_en);
    }

    cpgsave();
    cpgpage();
    printf("Finished\n\n");
}

int main(void)
{
    evaluate_ln_integral();
    root_finding_using_false_position();
    molecular_vibration();
}
