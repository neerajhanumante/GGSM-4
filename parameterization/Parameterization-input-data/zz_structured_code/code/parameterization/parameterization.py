'''
Geographical distribution related functions are populated here.
'''
from zz_structured_code.code.config.config_imports import *
from zz_structured_code.code.config.config_project_path import main_project_directory
from zz_structured_code.code.local_functions.local_functions import remove_border, f_read_df

def f_absolute_projections(df, df_country_codes, merge_column):
    
    # Merging the GDP and country-continent dfs on Country code
    df_var_countries_continents = pd.merge(left=df, right=df_country_codes, on=merge_column)
    # Aggregation
    df_var_countries_continents = df_var_countries_continents.groupby('Continent').sum()
    df_var_countries_continents.loc['World', :] = df_var_countries_continents.sum(axis=0)
    # Dropping unnecessary information
    try:
        df_var_countries_continents.drop('UN Code', inplace=True, axis=1)
        pass
    except:
        pass
    try:
        df_var_countries_continents.drop('Country code', inplace=True, axis=1)
        pass
    except:
        pass
    return df_var_countries_continents

    
def f_df_continents():
    '''
    UN code of a country and its location in a continent are obtained from the following source
    https://statisticstimes.com/geography/countries-by-continents.php accessed on 30 Sep 2021
    '''
    # Path creation
    directory_path = os.path.join(main_project_directory,
                                    'data',
                                    'input',)
    file_name = 'Country-continent-UN.csv'

    # Reading file
    df = pd.read_csv(join(directory_path, file_name), sep='\t', header=None)

    # Renaming columns
    l_col_headers = ['No', 'Country', 'Country code', 'UN Code', 'r1','r2', 'Continent']
    dct = {x:y for x, y in zip(np.arange(len(l_col_headers)), l_col_headers)}
    df.rename(columns=dct, inplace=True)

    # Selecting useful columns, sorting and resetting index
    df = df[['UN Code', 'Country', 'Continent', 'Country code']]
    df = df.sort_values(by=['Continent', 'UN Code'])
    df.reset_index(drop=True, inplace=True)

    return df


def f_extrapolate_backward_forward_df(df,
                                      extrapolate_upto_year=2120, 
                                      extrapolate_from_year=1950,
                                      data_year_range=[1940, 2060]):
    
    df_local = df.copy()

    def f_extend_df(from_year=1950, to_year=2120, df_local=df_local):
        complete_dataframe = df_local.copy()
        complete_dataframe.loc[:, np.arange(from_year, to_year+1)] = 0
        complete_dataframe.combine_first(df_local)
        return complete_dataframe

    def f_country_regressors():    
        condition_1 = [x for x in df_local.T.index 
                       if int(x) in np.arange(data_year_range[0], 
                                              data_year_range[1])]
        X = np.array(condition_1).reshape(-1, 1)
        l_regressors = []
    #     Obtaining regressors for each country
        for column_counter, column in enumerate(list(df_local.T.columns)):
#             print(column_counter, end=' ')
            y = np.array(df_local.loc[column, 
                                      np.array(condition_1)].T)
            regressor = LinearRegression()
            regressor.fit(X, y)
            l_regressors.append(regressor)
        return l_regressors
    
    
    def f_predict(prediction_horizon):
        prediction_horizon = prediction_horizon.reshape(-1, 1)
        # Prediction
        predicted_values = regressor.predict(prediction_horizon).tolist()
        for year_loc_counter, zipped_set in enumerate(zip(list(prediction_horizon), predicted_values)):
            year, value = zipped_set
            if value >= 0:
                last_positive_value = value
            if value <= 0:
            # What to do when projection becomes zero or negative?
            # There are two options: either 0 or last positive value
                if year_loc_counter > 0:
                    value = last_positive_value
                    value = 0
                else:
                    value = 0
                    last_positive_value = 0
            df_local.loc[column, year[0]] = value

    # Obtaining regressors 
    # Regressors are obtained for the df_local (before extending)
    l_regressors = f_country_regressors()
    
    # Extending the dataframe
    year_min = min(df_local.T.index)
    year_max = max(df_local.T.index)
    df_local = f_extend_df(extrapolate_from_year, extrapolate_upto_year, df_local)
    
    for column_counter, zipped_data in enumerate(zip(df_local.T.columns, l_regressors)):
        column, regressor = zipped_data
        # Creating new arrays for prediction
        prediction_horizon = np.arange(extrapolate_from_year,
                                       extrapolate_upto_year+1)
        f_predict(prediction_horizon)
#         prediction_horizon = np.arange(year_max+1,
#                                        extrapolate_upto_year+1)
#         f_predict(prediction_horizon)
    df_local.sort_index(axis=1, inplace=True)
    return df_local
