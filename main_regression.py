
from pathlib import Path
import numpy as np
import numpy.typing as npt

from curve_fitting import regression, plot

# Path to the Datasets directory
DATASETS_DIR = Path(__file__).parent / 'curve_fitting/data'

def f_b0(xi: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    return xi
def f_b(xi: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    return xi**2

def f(a: npt.NDArray[np.float64], xi: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    return a[0]*( 1.0 - np.exp(-a[1]*xi) )

def main():

    print('\nRegression')

    print('Linear Regression')
    file_path = DATASETS_DIR / 'linear_regression_data1.txt'
    data = np.loadtxt(file_path, delimiter=None, skiprows=0, usecols=(0,1))
    xi = data[:,0]
    yi = data[:,1]
    a = regression.linear_regression(xi, yi)
    print('a = ', a)
    xi_mat = xi[:, np.newaxis]
    y_model = regression.linear_model_predictions(xi_mat, a)
    print('y_model = ',y_model)
    n, m = xi.shape[0], 1
    std_error, r2_coef, ssr = regression.statistical_quantities(n, m, yi, y_model)
    print(f"Error: {std_error}, r^2: {r2_coef}, r: {r2_coef**0.5}, SSR: {ssr}")
    # plot.plot_experiment_model(xi, yi, y_model)

    print('\nPolynomial Regression')
    file_path = DATASETS_DIR / 'linear_regression_data2.txt'
    data = np.loadtxt(file_path, delimiter=None, skiprows=0, usecols=(0,1))
    xi = data[:,0]
    yi = data[:,1]
    a = regression.polynomial_regression(xi, yi, m=2)
    print('a = ', a)
    # Or
    xi_mat = regression.transform_data(xi, f_basis=[f_b0, f_b])
    a = regression.multi_linear_regression(xi_mat, yi)
    print('a = ', a)
    xi_mat = regression.transform_data(xi, f_basis=[f_b0, f_b])
    y_model = regression.linear_model_predictions(xi_mat, a)
    print('y_model = ',y_model)
    n, m = xi.shape[0], 1
    std_error, r2_coef, ssr = regression.statistical_quantities(n, m, yi, y_model)
    print(f"Error: {std_error}, r^2: {r2_coef}, r: {r2_coef**0.5}, SSR: {ssr}")
    # plot.plot_experiment_model(xi, yi, y_model)

    print('\nMultiple linear Regression')
    file_path = DATASETS_DIR / 'multi_linear_regression_data1.txt'
    data = np.loadtxt(file_path, delimiter=None, skiprows=0, usecols=(0,1,2))
    xi = data[:,:-1]
    yi = data[:,-1]
    a = regression.multi_linear_regression(xi, yi)
    print('a = ', a)
    y_model = regression.linear_model_predictions(xi, a)
    print('y_model = ',y_model)
    n, m = xi.shape[0], xi.shape[1]
    std_error, r2_coef, ssr = regression.statistical_quantities(n, m, yi, y_model)
    print(f"Error: {std_error}, r^2: {r2_coef}, r: {r2_coef**0.5}, SSR: {ssr}")
    # plot.plot_experiment_model(xi, yi, y_model)

    print('\nNon-linear Regression')
    file_path = DATASETS_DIR / 'non_linear_regression_data1.txt'
    data = np.loadtxt(file_path, delimiter=None, skiprows=0, usecols=(0,1))
    xi = data[:,0]
    yi = data[:,-1]
    a0 = np.array([1.0, 1.0])
    a = regression.non_linear_regression(f, xi, yi, a0, output=True)
    print(a)
    y_model = f(a, xi)
    n, m = xi.shape[0], 1
    std_error, r2_coef, ssr = regression.statistical_quantities(n, m, yi, y_model)
    print(f"Error: {std_error}, r^2: {r2_coef}, r: {r2_coef**0.5}, SSR: {ssr}")
    plot.plot_experiment_model(xi, yi, y_model)

if __name__ == '__main__':
    main()
