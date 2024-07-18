import traceback
import random
import sys
from datetime import datetime
import logging
import streamlit as st
import pandas as pd

from primertool.logger import init_logger

logger = init_logger(level=logging.ERROR, save_to="primertool.log")

widget_id = iter(random.sample(range(1, sys.maxsize), 1000))


def kuerzel_check(input_kuerzel: str) -> str:
    """Checks if the input initials are not empty and displays a warning toast if they are.

    Args:
        input_kuerzel (str): The input initials.

    Returns:
        str: The input initials.
    """
    if not input_kuerzel:
        st.toast(":orange[Don't forget to enter your initials]")
        input_kuerzel = None
    return input_kuerzel


def generate_primers(generator_function: callable, *args, **kwargs) -> pd.DataFrame or Exception:
    """Generates primers using the given generator function and its arguments and keyword arguments.

    Args:
        generator_function (callable): The primer generator function.
        *args: The arguments for the primer generator function.
        **kwargs: The keyword arguments for the primer generator function.

    Raises:
        Exception: If the primer generation fails.

    Returns:
        pd.DataFrame or Exception: The order table of the generated primers or an exception if the primer generation fails.
    """
    with st.spinner('Generating primers...'):
        try:
            return generator_function(*args, **kwargs).ordertable
        except Exception as e:
            logger.error(f"{traceback.format_exc()}")
            st.error(e)
            with st.expander("Advanced error information"):
                st.exception(e)


def feedback(message: str, start: datetime, df_ordertable: pd.DataFrame = None) -> None:
    """Displays a success message, the order table and the run time of the primer generation process.

    Args:
        message (str): The success message to display.
        start (datetime): The start time of the primer generation process.
        df_ordertable (pd.DataFrame, optional): The order table to display. Defaults to None.

    Returns:
        None
    """
    if df_ordertable is not None:
        st.success(message)
        st.dataframe(df_ordertable)
        runtime = datetime.now() - start
        st.info(f'Run time: {runtime.total_seconds():.2f} seconds')
