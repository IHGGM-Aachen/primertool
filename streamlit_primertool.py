from datetime import datetime
import streamlit as st

from primertool import primertool as pt
import streamlit_custom_components as cmp

widget_id = (hash(widget_id) for widget_id in range(1, 100_00))

st.markdown(f"This tool is built to facilitate the process of finding primers for PCR/Sanger sequencing. "
            f"Primers can be found for a given variant, exon, entire gene or genomic position.")
# ----------------------------------------------------------------------------------------------------------------------
#   EXPANDER
# ----------------------------------------------------------------------------------------------------------------------
with st.expander("How to use this tool"):
    st.markdown("To get started, please")
    st.markdown("- enter your initials (Kürzel); Default is None")
    st.markdown(f"- select the genome assembly you are working with (default is hg38)")
    st.markdown(f"- select the primer generation method that corresponds to the type of primer you are looking for and "
                f"enter the required information")
    st.markdown(f'- hit the "Generate Primers 🧬" button and wait until the order table is generated (depending on your '
                f'input this may take some time, but usually less than a minute)')
    st.markdown(f"Note: When running the tool for the first time, it may take several minutes to download the selected "
                f"reference genome. Subsequent runs will be faster.")

# ----------------------------------------------------------------------------------------------------------------------
#   INPUTS (Initials & Genome Assembly)
# ----------------------------------------------------------------------------------------------------------------------
col_kuerzel, col_genome_assembly = st.columns(2)
with col_kuerzel:
    kuerzel = st.text_input(label='Initials', placeholder='AB', help='Enter your initials (e.g. "XY")',
                            key=next(widget_id))
with col_genome_assembly:
    genome_assembly = st.selectbox('Genome assembly', ('hg38', 'hg19'), index=0, help='Select genome assembly',
                                   key=next(widget_id))

# ----------------------------------------------------------------------------------------------------------------------
#   TAB MENU
# ----------------------------------------------------------------------------------------------------------------------
tab_variant, tab_exon, tab_gene, tab_genomic_pos = st.tabs(["Variant", "Exon", "Gene", "Genomic Position"])

# ----------------------------------------------------------------------------------------------------------------------
# ----- VARIANT PRIMER GENERATOR
# ----------------------------------------------------------------------------------------------------------------------
with tab_variant, st.form('variant_form', border=False):
    st.subheader("Generate primers for a variant")
    hgvs_variant = st.text_input(label='Variant (e.g. "NM_000410.3\:c.845G>A")', placeholder='NM_000410.3:c.845G>A',
                                 help='Enter variant in [HGVS notation](https://hgvs-nomenclature.org/stable/) (e.g. "NM_000410.3\:c.845G>A")',
                                 key=next(widget_id))
    submitted = st.form_submit_button('Generate Primers 🧬', type='primary')

    if submitted and genome_assembly and hgvs_variant:
        st.divider()
        kuerzel = cmp.kuerzel_check(kuerzel)
        start_time = datetime.now()
        df_ordertable = cmp.generate_primers(pt.VariantPrimerGenerator, hgvs_variant, genome_assembly, kuerzel=kuerzel)
        hgvs_variant = hgvs_variant.replace(':', '\:')
        cmp.feedback(f"Primers for variant {hgvs_variant} in genome assembly {genome_assembly} are:",
                     start_time, df_ordertable)
    elif submitted:
        st.warning("Please enter a variant and select a genome assembly")
# ----------------------------------------------------------------------------------------------------------------------
# ----- EXON PRIMER GENERATOR
# ----------------------------------------------------------------------------------------------------------------------
with tab_exon, st.form('exon_form', border=False):
    st.subheader("Generate primers for an exon")
    nm_number = st.text_input(label='NM number', placeholder='NM_000451', help='Enter NM number (e.g. "NM_000451")',
                              key=next(widget_id))
    exon_number = st.number_input(label='Exon number', min_value=0, step=1, help='Enter exon number',
                                  key=next(widget_id))
    submitted = st.form_submit_button('Generate Primers 🧬', type='primary')

    if submitted and genome_assembly and nm_number and exon_number:
        st.divider()
        kuerzel = cmp.kuerzel_check(kuerzel)

        start_time = datetime.now()
        df_ordertable = cmp.generate_primers(pt.ExonPrimerGenerator, nm_number, exon_number, genome_assembly,
                                             kuerzel=kuerzel)
        cmp.feedback(f"Primers for exon {exon_number} in gene {nm_number} in genome assembly {genome_assembly} "
                     f"are:", start_time, df_ordertable)
    elif submitted:
        st.warning("Please enter NM number, exon number and select a genome assembly")
# ----------------------------------------------------------------------------------------------------------------------
# ----- GENE PRIMER GENERATOR
# ----------------------------------------------------------------------------------------------------------------------
with tab_gene, st.form('gene_form', border=False):
    st.subheader("Generate primers for a gene")
    nm_number = st.text_input(label='NM number', placeholder='NM_000451', help='Enter NM number (e.g. "NM_000451")',
                              key=next(widget_id))
    submitted = st.form_submit_button('Generate Primers 🧬', type='primary')

    if submitted and genome_assembly and nm_number:
        st.divider()
        kuerzel = cmp.kuerzel_check(kuerzel)

        start_time = datetime.now()
        df_ordertable = cmp.generate_primers(pt.GenePrimerGenerator, nm_number, genome_assembly, kuerzel=kuerzel)
        cmp.feedback(f"Primers for gene {nm_number} in genome assembly {genome_assembly} are:",
                     start_time, df_ordertable)
    elif submitted:
        st.warning("Please enter NM number and select a genome assembly")
# ----------------------------------------------------------------------------------------------------------------------
# ----- GENOMIC POSITION PRIMER GENERATOR
# ----------------------------------------------------------------------------------------------------------------------
with tab_genomic_pos, st.form('genomic_pos_form', border=False):
    st.subheader("Generate primers for a genomic position")

    chromosome = st.text_input(label='Chromosome', placeholder='chr19', help='Enter chromosome (e.g. "chr1")',
                               key=next(widget_id))
    start_position = st.number_input(label='Start', min_value=0, step=1, help='Enter start position', key=next(widget_id))
    end_position = st.number_input(label='End', min_value=0, step=1, help='Enter end position', key=next(widget_id))

    # optional input container
    container = st.container(border=True)
    container.markdown(":grey[Optional input: only appears on order table, not relevant for computation]")
    col_gene, col_nm_number = container.columns(2)
    with col_gene:
        gene = st.text_input(label='Gene (optional)', placeholder='BRCA1', help='Enter gene (e.g. "BRCA1")',
                             key=next(widget_id))
    with col_nm_number:
        nm_number = st.text_input(label='NM number (optional)', placeholder='NM_000451',
                                  help='Enter NM number (e.g. "NM_000451")', key=next(widget_id))

    submitted = st.form_submit_button('Generate Primers 🧬', type='primary')

    if submitted and genome_assembly and chromosome and start_position and end_position:
        st.divider()
        kuerzel = cmp.kuerzel_check(kuerzel)

        start_time = datetime.now()
        df_ordertable = cmp.generate_primers(pt.GenomicPositionPrimerGenerator, chromosome, start_position,
                                             end_position, genome_assembly, kuerzel=kuerzel)
        if df_ordertable is not None and not df_ordertable.empty and gene:
            df_ordertable['gene'] = gene
        if df_ordertable is not None and not df_ordertable.empty and nm_number:
            df_ordertable['nm_number'] = nm_number

        cmp.feedback(f"Primers for genomic position {chromosome}:{start_position}-{end_position} in genome "
                     f"assembly {genome_assembly} are:", start_time, df_ordertable)
    elif submitted:
        st.warning("Please enter chromosome, start position, end position and select a genome assembly")
