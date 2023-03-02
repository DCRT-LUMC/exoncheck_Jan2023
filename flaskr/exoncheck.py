from flask import Blueprint
from flask import flash
from flask import g
from flask import redirect
from flask import request
from flask import url_for
from werkzeug.exceptions import abort
from flaskr.auth import login_required
from flaskr.db import get_db
from flaskr.info import *
from flask import render_template
from flask import session

bp = Blueprint("exoncheck", __name__)

@bp.route("/")
def welcome():
    """Show the welcome page as default/root."""
    return render_template("exoncheck/welcome.html")

@bp.route("/information")
def information():
    """Show additional information about eligibility metrics."""
    return render_template("exoncheck/information.html")

@bp.route("/output")
@login_required
def output():
    """Show all uploaded variants of the current user, most recent first"""
    db = get_db()

    variants = db.execute(
        "SELECT v.id, input_variant, created, author_id, username, "

        # Show hg info
        "hg,"

        # Show gene info
        "gene_symbol,"
        "ENSG_gene_id,"
        "elig,"

        # Show protein info
        "prot_name,"
        "short_name,"

        # Show variant info
        "NC_variant,"
        "strand,"
        "hg_variant,"
        "MANE_select_NM_variant,"
        "MANE_select_ENST_variant,"
        "consequence_variant,"

        # Show exon info
        "exon_number,"
        "total_exons,"
        "exon_number_interpretation,"
        "coding_exons,"
        "NC_exon,"
        "exon_length,"
        "nearest_splice_distant,"
        "nearest_end,"
        "total_protein_length,"
        "percentage_length,"
        "length_condition,"
        "frame,"
        "splice_dist_interpretation,"
        "MANE_select_NM_exon,"
        "r_exon_skip,"
#        "MANE_select_ENST_exon,"
        "consequence_skipping,"

        # Show domain info
        "uniprot_link,"
        "domain_info,"

        # Show expression info
        "gtex_link,"
        "expression_brain,"
        "expression_fibroblasts,"
        "expression_tibial_nerve,"
        "expression_blood,"
        "expression_transformed_lymphocytes,"
        "expression_periphery_retina,"
        "expression_center_retina,"

        # Show LOVD info
        "lovd_output,"

        # Show OMIM id
        "omim_id,"

        # Show links
        "omim_link,"
        "gnomAD_link,"
        "decipher_link,"
        "clinvar_link"

        # Show uploads of the current user sorted by recent uploads first
        " FROM variant v JOIN user u ON v.author_id = u.id"
        " WHERE u.id = ? ORDER BY created DESC", (g.user["id"],),
    ).fetchall()
    return render_template("exoncheck/output.html", variants=variants)

@bp.route("/gene-output")
@login_required
def gene_output():
    """ Show all selected exons of the input gene of the current user, most recent first"""
    db = get_db()

    exons = db.execute(
        "SELECT e.id, input_gene, created, author_id, "
        "NM_id, NC_id, elig, "
        "exon_number, total_exons, exon_number_interpretation, "
        "coding_exons, NC_exon, exon_length, "
        "total_protein_length, percentage_length, length_condition, "
        "frame, MANE_select_NM_exon, r_exon_skip, consequence_skipping "

        # Show uploads of the current user sorted by recent uploads first
        " FROM exon e JOIN user u ON e.author_id = u.id"
        " WHERE u.id = ? ORDER BY created DESC", (g.user["id"],),
        ).fetchall()
    return render_template("exoncheck/gene_output.html", exons=exons)

def get_variant(id, check_author=True):
    """
    Get an uploaded variant and its author by id
    Checks that the id exists
    Checks if current user is the author

    :param id: id of post to get
    :param check_author: require the current user to be the author
    :return: the post with author information
    :raise 404: if a post with the given id doesn't exist
    :raise 403: if the current user isn't the author
    """
    variant = (
        get_db()
        .execute(
            "SELECT v.id, input_variant, created, author_id, username"
            " FROM variant v JOIN user u ON v.author_id = u.id"
            " WHERE v.id = ?",
            (id,),
        )
        .fetchone()
    )

    if variant is None:
        abort(404, f"Variant id {id} doesn't exist.")

    if check_author and variant["author_id"] != g.user["id"]:
        abort(403)
    return variant

def get_exon(id, check_author=True):
    """
    Get a selected exon and its author by id
    Checks that the id exists
    Checks if current user is the author

    :param id: id of exon to get
    :param check_author: require the current user to be the author
    :return: the exon with author information
    :raise 404: if a post with the given id doesn't exist
    :raise 403: if the current user isn't the author
    """
    exon = (
        get_db()
        .execute(
            "SELECT e.id, exon_number, created, author_id, username"
            " FROM exon e JOIN user u ON e.author_id = u.id"
            " WHERE e.id = ?",
            (id,),
        )
        .fetchone()
    )

    if exon is None:
        abort(404, f"Exon id {id} doesn't exist.")

    if check_author and exon["author_id"] != g.user["id"]:
        abort(403)
    return exon

@bp.route("/exons", methods=("GET", "POST"))
@login_required
def exons():
    input_gene = session.get("input_gene")
    NM_id, NC_id, in_frame_exons, out_of_frame_exons = gene_to_exons(input_gene)
    elig = check_gene_eligibility(input_gene)

    if request.method == "POST":
        exon_number = request.form.get("exon_number")

        total_exons,\
        total_protein_length,\
        coding_exons,\
        NC_exon,\
        coding_exon_length,\
        exon_number_interpretation,\
        exon_length,\
        frame,\
        percentage_length,\
        length_condition,\
        MANE_select_NM_exon,\
        consequence_skipping,\
        r_exon_skip = exon_skip(NM_id, NC_id, exon_number)

        db = get_db()
        db.execute(
        "INSERT INTO exon "
        # user
        "(input_gene,"
        "author_id,"

        # gene
        "NM_id,"
        "NC_id,"
        "elig,"

        # exon
        "exon_number,"
        "total_exons,"
        "exon_number_interpretation,"
        "coding_exons,"
        "NC_exon,"
        "exon_length,"
        "total_protein_length,"
        "percentage_length,"
        "length_condition,"
        "frame,"
        "MANE_select_NM_exon,"
        "r_exon_skip,"
        "consequence_skipping)"
        "VALUES (?, ?, ?, ?, ?,"
        "?, ?, ?, ?, ?, ?, ?, ?,"
        "?, ?, ?, ?, ?)",

        (# user
        input_gene,
        g.user["id"],

        # gene
        NM_id,
        NC_id,
        elig,

        # exon
        exon_number,
        total_exons,
        exon_number_interpretation,
        coding_exons,
        NC_exon,
        exon_length,
        total_protein_length,
        percentage_length,
        length_condition,
        frame,
        MANE_select_NM_exon,
        r_exon_skip,
        consequence_skipping))

        db.commit()
        return redirect(url_for("exoncheck.gene_output"))
    return render_template("exoncheck/exons.html", NM_id = NM_id, NC_id = NC_id, in_frame_exons = in_frame_exons, out_of_frame_exons = out_of_frame_exons)

@bp.route("/create", methods=("GET", "POST"))
@login_required
def create():
    """Upload new variant by the current user."""
    if request.method == "POST":
        input_variant = request.form.get("input_variant")
        input_gene = request.form.get("input_gene")
        hg = request.form.get("hg")
        syntax_message = check_for_hgvs_format(input_variant)

        if input_gene is not None:
            session["input_gene"] = input_gene
            return redirect(url_for("exoncheck.exons", input_gene = input_gene))

        if input_variant is not None:
            if input_variant.startswith('NP'):
                message, NM_variants = get_NM_for_NP(input_variant)
                if message != "":
                    flash(message)
                else:
        #                db = get_db()
        #                db.execute("INSERT INTO post (NM_variants) VALUES (?);", (NM_variants))
        #                db.commit()
        #                NM_variants = json.loads(post['NM_variants'])
                    return render_template("exoncheck/variants.html", NM_variants = NM_variants)
            else:
                if syntax_message != '':
                    flash(syntax_message)

                else:
                    # Check if input is in correct HGVS format using VariantValidator, syntax only
                    # Check if nucleotides and positions of input correspond to the submitted transcript
                    match_message = check_for_match_variant_and_transcript(input_variant)

                    if match_message != '':
                        flash(match_message)

                    else:
                        # Get MANE select in NM format
                        MANE_select_NM_variant, \
                        MANE_select_ENST_variant = get_MANE_select_identifiers(input_variant)

                        if hg == 'hg19/GRCh37':
                            # Exploit variant validator
                            NC_variant, \
                            hg_variant, \
                            ENSG_gene_id, \
                            omim_id, \
                            gene_symbol, \
                            consequence_variant, \
                            exon_number, \
                            total_exons, \
                            exon_number_interpretation, \
                            coding_exons, \
                            NC_exon, \
                            exon_length, \
                            nearest_splice_distant, \
                            nearest_end, \
                            total_protein_length, \
                            percentage_length, \
                            length_condition, \
                            frame, \
                            splice_dist_interpretation, \
                            consequence_skipping, \
                            r_exon_skip, \
                            MANE_select_NM_exon = exploit_variant_validator_hg19(input_variant)

                            # Get gene eligibility
                            elig = check_gene_eligibility(gene_symbol)

                            # Get LOVD information
                            lovd_output = get_lovd_info_hg19(hg_variant, NC_variant)

                            # Get strand
                            strand = get_strand(ENSG_gene_id)

                            # Get links
                            gtex_link = 'https://gtexportal.org/home/gene/' + ENSG_gene_id
                            omim_link = 'https://www.omim.org/entry/' + omim_id
                            gnomAD_link = 'https://gnomad.broadinstitute.org/gene/' + ENSG_gene_id + '?dataset=gnomad_r2_1'
                            decipher_link = 'https://www.deciphergenomics.org/sequence-variant/' + hg_variant
                            clinvar_link = 'https://www.ncbi.nlm.nih.gov/clinvar/?term=' + MANE_select_NM_variant

                            uniprot_link, \
                            domain_info, \
                            prot_name, \
                            short_name = get_uniprot_info(ENSG_gene_id, r_exon_skip)

                            # Get expressions
                            expression_brain, expression_fibroblasts, expression_tibial_nerve, expression_blood, expression_transformed_lymphocytes = get_gene_expression(ENSG_gene_id, MANE_select_ENST_variant)
                            expression_periphery_retina, expression_center_retina = get_eye_expression(ENSG_gene_id)

                        else:
                            # Exploit variant validator
                            NC_variant, \
                            hg_variant, \
                            ENSG_gene_id, \
                            omim_id, \
                            gene_symbol, \
                            consequence_variant, \
                            exon_number, \
                            total_exons, \
                            exon_number_interpretation, \
                            coding_exons, \
                            NC_exon, \
                            exon_length, \
                            nearest_splice_distant, \
                            nearest_end, \
                            total_protein_length, \
                            percentage_length, \
                            length_condition, \
                            frame, \
                            splice_dist_interpretation, \
                            consequence_skipping, \
                            r_exon_skip, \
                            MANE_select_NM_exon = exploit_variant_validator(MANE_select_NM_variant)

                            # Get gene eligibility
                            elig = check_gene_eligibility(gene_symbol)

                            # Get LOVD information
                            lovd_output = get_lovd_info(hg_variant, NC_variant)

                            # Get strand
                            strand = get_strand(ENSG_gene_id)

                            # Get links
                            gtex_link = 'https://gtexportal.org/home/gene/' + ENSG_gene_id
                            omim_link = 'https://www.omim.org/entry/' + omim_id
                            gnomAD_link = 'https://gnomad.broadinstitute.org/gene/' + ENSG_gene_id + '?dataset=gnomad_r3'
                            decipher_link = 'https://www.deciphergenomics.org/sequence-variant/' + hg_variant
                            clinvar_link = 'https://www.ncbi.nlm.nih.gov/clinvar/?term=' + MANE_select_NM_variant

                            uniprot_link, \
                            domain_info, \
                            prot_name, \
                            short_name = get_uniprot_info(ENSG_gene_id, r_exon_skip)

                            # Get expressions
                            expression_brain, expression_fibroblasts, expression_tibial_nerve, expression_blood, expression_transformed_lymphocytes = get_gene_expression(ENSG_gene_id, MANE_select_ENST_variant)
                            expression_periphery_retina, expression_center_retina = get_eye_expression(ENSG_gene_id)


                        db = get_db()
                        db.execute(
                            "INSERT INTO variant "
                            # user
                            "(input_variant,"
                            "author_id,"

                            # hg-version
                            "hg,"

                            # gene
                            "gene_symbol,"
                            "ENSG_gene_id,"
                            "elig,"

                            # protein
                            "prot_name,"
                            "short_name,"

                            # variant
                            "NC_variant,"
                            "strand,"
                            "hg_variant,"
                            "MANE_select_NM_variant,"
                            "MANE_select_ENST_variant,"
                            "consequence_variant,"

                            # exon
                            "exon_number,"
                            "total_exons,"
                            "exon_number_interpretation,"
                            "coding_exons,"
                            "NC_exon,"
                            "exon_length,"
                            "nearest_splice_distant,"
                            "nearest_end,"
                            "total_protein_length,"
                            "percentage_length,"
                            "length_condition,"
                            "frame,"
                            "splice_dist_interpretation,"
                            "MANE_select_NM_exon,"
                            "r_exon_skip,"
    #                        "MANE_select_ENST_exon,"
                            "consequence_skipping,"

                            # domain
                            "uniprot_link,"
                            "domain_info,"

                            # expression
                            "gtex_link,"
                            "expression_brain,"
                            "expression_fibroblasts,"
                            "expression_tibial_nerve,"
                            "expression_blood,"
                            "expression_transformed_lymphocytes,"
                            "expression_periphery_retina,"
                            "expression_center_retina,"

                            # lovd
                            "lovd_output,"

                            # identifiers
                            "omim_id,"

                            # links
                            "omim_link,"
                            "gnomAD_link,"
                            "decipher_link,"
                            "clinvar_link)"
                            "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?,"
                            "?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?,"
                            "?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?,"
                            "?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",

                            (# user
                            input_variant,
                            g.user["id"],

                            # hg-version
                            hg,

                            # gene
                            gene_symbol,
                            ENSG_gene_id,
                            elig,

                            # protein
                            prot_name,
                            short_name,

                            # variant
                            NC_variant,
                            strand,
                            hg_variant,
                            MANE_select_NM_variant,
                            MANE_select_ENST_variant,
                            consequence_variant,

                            # exon
                            exon_number,
                            total_exons,
                            exon_number_interpretation,
                            coding_exons,
                            NC_exon,
                            exon_length,
                            nearest_splice_distant,
                            nearest_end,
                            total_protein_length,
                            percentage_length,
                            length_condition,
                            frame,
                            splice_dist_interpretation,
                            MANE_select_NM_exon,
                            r_exon_skip,
    #                        MANE_select_ENST_exon,
                            consequence_skipping,

                            # domain
                            uniprot_link,
                            domain_info,

                            # expression
                            gtex_link,
                            expression_brain,
                            expression_fibroblasts,
                            expression_tibial_nerve,
                            expression_blood,
                            expression_transformed_lymphocytes,
                            expression_periphery_retina,
                            expression_center_retina,

                            # lovd
                            lovd_output,

                            # identifiers
                            omim_id,

                            # links
                            omim_link,
                            gnomAD_link,
                            decipher_link,
                            clinvar_link))

                        db.commit()
                        return redirect(url_for("exoncheck.output"))
    return render_template("exoncheck/create.html")

@bp.route("/<int:id>/update", methods=("GET", "POST"))
@login_required
def update(id):
    """Update a post if the current user is the author."""
    variant = get_variant(id)

    if request.method == "POST":
        input_variant = request.form["input_variant"]
        body = request.form["body"]
        error = None

        if not input_variant:
            error = "An input is required."

        if error is not None:
            flash(error)
        else:
            db = get_db()
            db.execute(
                "UPDATE variant SET input_variant = ?, body = ? WHERE id = ?", (input_variant, body, id)
            )
            db.commit()
            return redirect(url_for("exoncheck.output"))

    return render_template("exoncheck/update.html", variant=variant)

@bp.route("/<int:id>/delete", methods=("POST",))
@login_required
def delete(id):
    """Delete a post.
    Ensures that the post exists and that the logged in user is the
    author of the post.
    """
    get_variant(id)
    db = get_db()
    db.execute("DELETE FROM variant WHERE id = ?", (id,))
    db.commit()
    return redirect(url_for("exoncheck.output"))

@bp.route("/<int:id>/updateExon", methods=("GET", "POST"))
@login_required
def updateExon(id):
    """Update an exon if the current user is the author."""
    exon = get_exon(id)

    if request.method == "POST":
        exon_number = request.form["exon"]
        body = request.form["body"]
        error = None

        if not exon_number:
            error = "An input is required."

        if error is not None:
            flash(error)
        else:
            db = get_db()
            db.execute(
                "UPDATE exon SET exon_number = ?, body = ? WHERE id = ?", (exon_number, body, id)
            )
            db.commit()
            return redirect(url_for("exoncheck/gene_output.html"))

    return render_template("exoncheck/update_exon.html", exon=exon)


@bp.route("/<int:id>/deleteExon", methods=("POST",))
@login_required
def deleteExon(id):
    """Delete a post.
    Ensures that the post exists and that the logged in user is the
    author of the post.
    """
    get_exon(id)
    db = get_db()
    db.execute("DELETE FROM exon WHERE id = ?", (id,))
    db.commit()
    return redirect(url_for("exoncheck.gene_output"))
