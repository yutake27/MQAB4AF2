"""
Function for MQA performance evaluation
"""

from typing import List, Union
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from sklearn.metrics import mean_absolute_error, mean_squared_error

pd.options.display.float_format = '{:.3f}'.format
plt.rcParams["figure.dpi"] = 300
sns.set(style='darkgrid')
warnings.simplefilter('ignore', UserWarning)


def get_loss(target_df: pd.DataFrame, label_name: str, column_name: str, how='mean') -> Union[float, None]:
    """
        target_df: DataFrame of a target.
        label_name: Label of model quality. e.g) GDT_TS, GDT_HA
        column_name: A method name.
        how: {'mean', 'best', 'random'}, default: 'mean'
    """
    assert how in ['mean', 'best', 'random']
    try:
        column_max_value = target_df[column_name].max()
        if how == 'mean':
            column_max_value_df = target_df[target_df[column_name] == column_max_value]
            loss = target_df[label_name].max() - column_max_value_df[label_name].mean()
        elif how == 'best':
            column_max_value_df = target_df[target_df[column_name] == column_max_value]
            loss = target_df[label_name].max() - column_max_value_df[label_name].max()
        elif how == 'random':
            loss = target_df[label_name].max() - target_df[label_name][target_df[column_name].idxmax()]
    except (TypeError, KeyError):
        loss = None
    return loss


def get_whole_loss(df: pd.DataFrame, label_name: str, column_list: list, how: str = 'mean') -> pd.DataFrame:
    loss_list = [get_loss(df, label_name, column, how=how) for column in column_list]
    return pd.DataFrame({label_name + ' Loss': loss_list}, index=column_list)


def get_mae(df: pd.DataFrame, label_name: str, column_name: str) -> Union[float, None]:
    df = df.dropna(subset=[column_name])
    if df[column_name].max() > 1:
        return None
    try:
        mae = mean_absolute_error(df[label_name], df[column_name])
    except (ValueError, KeyError):
        mae = None
    return mae


def get_whole_mae(df: pd.DataFrame, label_name: str, column_list: list) -> pd.DataFrame:
    mae_list = [get_mae(df, label_name, column) for column in column_list]
    return pd.DataFrame({label_name + ' MAE': mae_list}, index=column_list)


def get_mse(df: pd.DataFrame, label_name: str, column_name: str) -> Union[float, None]:
    df = df.dropna(subset=[column_name])
    try:
        mse = mean_squared_error(df[label_name], df[column_name])
    except (ValueError, KeyError):
        mse = None
    return mse


def get_whole_mse(df: pd.DataFrame, label_name: str, column_list: list) -> pd.DataFrame:
    mse_list = [get_mse(df, label_name, column) for column in column_list]
    return pd.DataFrame({label_name + ' MSE': mse_list}, index=column_list)


def eval_get_df(df: pd.DataFrame, columns: List, label_name: str = 'GDT_TS',
                threshold: float = 0, loss_how: str = 'mean') -> pd.DataFrame:
    """Evaluate MQA performance for each target in a DataFrame.

    Args:
        df (pd.DataFrame): DataFrame of MQA scores.
        columns (List): List of method names.
        label_name (str, optional): Label name. Defaults to 'GDT_TS'.
        threshold (float, optional): Min threshold to remove targets. Defaults to 0.
        loss_how (str, optional): {'mean', 'best', 'random'}. Defaults to 'mean'.

    Returns:
        pd.DataFrame: DataFrame of MQA performance for each target.
    """
    df = df.groupby('Target').filter(lambda x: x[label_name].max() >= threshold)
    if df[label_name].max() > 1:
        df[label_name] /= 100
    group = df.groupby('Target')
    pearson = group.corr()[label_name].loc[:, columns].rename(label_name + ' Pearson')
    spearman = group.corr(method='spearman')[label_name].loc[:, columns].rename(label_name + ' Spearman')
    loss = group.apply(lambda x: get_whole_loss(x, label_name, columns, how=loss_how)) * 100
    mae = group.apply(lambda x: get_whole_mae(x, label_name, columns))
    pef_df = pd.concat([pearson, spearman, loss, mae], axis=1).reset_index().rename(columns={'level_1': 'Method'})
    print(len(group))
    return pef_df


def eval(df: pd.DataFrame, columns: List, label_name: str = 'GDT_TS', threshold: float = 0, **kwargs) -> pd.DataFrame:
    """Evaluate MQA performance for all targets in a DataFrame.

    Args:
        df (pd.DataFrame): DataFrame of MQA scores.
        columns (List): List of method names.
        label_name (str, optional): Label names. Defaults to 'GDT_TS'.
        threshold (float, optional): Threshold to remove the target.
        Remove targets whose maximum label value is less than the threshold.　Defaults to 0.

    Returns:
        pd.DataFrame: DataFrame of MQA performance for the dataset.
    """
    pef_df = eval_get_df(df, columns=columns, label_name=label_name, threshold=threshold,
                         **kwargs).groupby('Method').mean().reset_index()
    order_dict = dict(zip(columns, range(len(columns))))
    pef_df['order'] = [order_dict[method] for method in pef_df['Method']]
    pef_df = pef_df.sort_values('order')
    pef_df = pef_df.set_index('Method')
    pef_df = pef_df.reset_index().drop('order', axis=1)
    return pef_df


def stat_test(df, base_method='pLDDT', methods=None, metrics=None):
    """Conduct statistical test

    Args:
        df (pd.DataFrame): DataFrame of MQA performance for each target
        base_method (str, optional): Name of the base method. Defaults to 'identity(%)'.
        methods (list, optional): List of the methods. Defaults to None.
        metrics (list, optional): List of the metrics. Defaults to None.

    Returns:
        pd.DataFrame: DataFrame of the result of the statistical test.
    """
    if methods is None:
        methods = ['DOPE', 'SOAP', 'ProQ3D', 'SBROD', 'DeepAccNet', 'DeepAccNet-Bert',
                   'P3CMQA', 'VoroCNN', 'pLDDT', 'pTMscore']
    if metrics is None:
        prefix = 'GDT_TS '
        metrics = [prefix + metric for metric in ['Pearson', 'Spearman', 'Loss', 'MAE']]

    df = df.sort_values('Target')
    p_dic_list = []
    metric_prefix = ''
    # for each method
    for method in methods:
        # skip the base method
        if method == base_method:
            p_dic = {'Method': base_method}
            for metric in metrics:
                p_dic[metric_prefix + metric] = None
            p_dic_list.append(p_dic)
            continue
        # for other methods
        p_dic = {}
        p_dic['Method'] = method
        for metric in metrics:
            base_metric_score = df.query('Method == @base_method')[metric]
            comp_metric_score = df.query('Method == @method')[metric]
            if comp_metric_score.isnull().all():
                metrics_p_value = None
            else:
                metrics_p_value = stats.wilcoxon(base_metric_score, comp_metric_score)[1]
            p_dic[metric_prefix + metric] = metrics_p_value
        p_dic_list.append(p_dic)
    p_value_df = pd.DataFrame(p_dic_list)
    return p_value_df
