import streamlit as st

# ----------------------------------------------------------------------------------------------------------------------
#   HEADER
# ----------------------------------------------------------------------------------------------------------------------
st.set_page_config(page_title="Primertool", page_icon="logo.png")  # Dna icons created by Freepik - Flaticon

# ----------------------------------------------------------------------------------------------------------------------
#   TITLE
# ----------------------------------------------------------------------------------------------------------------------
st.title("Primertool")

main_page = st.Page('streamlit_primertool.py', title='Main', icon=':material/genetics:')
help_page = st.Page('streamlit_help.py', title='Help', icon=':material/question_mark:')
pg = st.navigation([main_page, help_page])
pg.run()

# ----------------------------------------------------------------------------------------------------------------------
# st.balloons()
# st.snow()
