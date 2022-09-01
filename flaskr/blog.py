from flask import Blueprint
from flask import flash
from flask import g
from flask import redirect
from flask import render_template
from flask import request
from flask import url_for
from werkzeug.exceptions import abort

from flaskr.auth import login_required
from flaskr.db import get_db
from flaskr.info import *

bp = Blueprint("blog", __name__)

@bp.route("/welcome")
def welcome():
    return render_template("blog/welcome.html")


@bp.route("/output")
def output():
    """Use fetchall() for: Show all the posts, most recent first.
    But we will now just show the latest variant"""
    db = get_db()

    posts = db.execute(
        "SELECT p.id, title, created, author_id, username, "
        # gene
        "gene_symbol,"
        "ENSG_gene_id,"

        # variant
        "NC_variant,"
        "strand,"
        "hg38_variant,"
        "MANE_select_NM_variant,"
        "MANE_select_ENST_variant,"
        "consequence_variant,"

        # exon
        "exon_number,"
        "total_exons,"
        "NC_exon,"
        "exon_length,"
        "total_protein_length,"
        "percentage_length,"
        "frame,"
        "MANE_select_NM_exon,"
        "MANE_select_ENST_exon,"
        "consequence_skipping,"

        # domain
        "uniprot_link,"
        "domain_info,"

        # expression
        "gtex_link,"
        "expression_eye,"
        "expression_brain,"
        "expression_fibroblasts,"
        "expression_tibial_nerve,"
        "expression_blood,"
        "expression_transformed_lymphocytes,"

        # lovd
        "no_exact_lovd_matches,"
        "no_partial_lovd_matches,"
        "exact_lovd_match_link,"
        "lovd_link,"
        "partial_lovd_match1,"
        "partial_lovd_match2,"
        "partial_lovd_match3,"
        "partial_lovd_match4,"
        "partial_lovd_match5,"
        "partial_lovd_match6,"
        "partial_lovd_match7,"
        "partial_lovd_match8,"
        "partial_lovd_match9,"
        "partial_lovd_match10,"

        # identifiers
        "omim_id,"

        # links
        "omim_link,"
        "gnomAD_link,"
        "decipher_link,"
        "clinvar_link"

        " FROM post p JOIN user u ON p.author_id = u.id"
        " ORDER BY created DESC"
    ).fetchmany()
    return render_template("blog/output.html", posts=posts)


def get_post(id, check_author=True):
    """Get a post and its author by id.

    Checks that the id exists and optionally that the current user is
    the author.

    :param id: id of post to get
    :param check_author: require the current user to be the author
    :return: the post with author information
    :raise 404: if a post with the given id doesn't exist
    :raise 403: if the current user isn't the author
    """
    post = (
        get_db()
        .execute(
            "SELECT p.id, title, created, author_id, username"
            " FROM post p JOIN user u ON p.author_id = u.id"
            " WHERE p.id = ?",
            (id,),
        )
        .fetchone()
    )

    if post is None:
        abort(404, f"Post id {id} doesn't exist.")

    if check_author and post["author_id"] != g.user["id"]:
        abort(403)

    return post


@bp.route("/create", methods=("GET", "POST"))
@login_required
def create():
    """Create a new post for the current user."""
    if request.method == "POST":
        title = request.form["title"]
        error_title = None

        if not title:
            error_title = "Please upload your variant in HGVS format"

        if error_title is not None:
            flash(error_title)

        # Check if input is in correct HGVS format using VariantValidator
        else:
            error_hgvs_format = None
            if check_for_hgvs_format(title) != "":
                error_hgvs_format = f"Your variant ({title}) is not conform HGVS format. " \
                                    f"See http://varnomen.hgvs.org/recommendations/general/ for a quick overview. " \
                                    f"Example input: NM_001029896.2:c.749_751del. " \
                                    f"Please try again." \

                                    # Optionally: show error message
                                    # "(Feedback message =\t" + check_for_hgvs_format(title) + ")"

            if error_hgvs_format is not None:
                flash(error_hgvs_format)

            else:
                # Get MANE select in NM format
                MANE_select_NM_variant, \
                MANE_select_ENST_variant = get_MANE_select_identifiers(title)

                # Exploit variant validator
                NC_variant, \
                hg38_variant, \
                ENSG_gene_id, \
                omim_id, \
                gene_symbol, \
                consequence_variant, \
                exon_number, \
                total_exons, \
                NC_exon, \
                exon_length, \
                total_protein_length, \
                percentage_length, \
                frame, \
                consequence_skipping, \
                MANE_select_NM_exon = exploit_variant_validator(MANE_select_NM_variant)

                # Get LOVD information
                lovd_link, no_exact_lovd_matches, exact_lovd_match_link, no_partial_lovd_matches, partial_lovd_match1, partial_lovd_match2, partial_lovd_match3, partial_lovd_match4, \
                partial_lovd_match5, partial_lovd_match6, partial_lovd_match7, partial_lovd_match8, \
                partial_lovd_match9, partial_lovd_match10 = get_lovd_info(NC_variant, gene_symbol)

                # Get strand
                strand = get_strand(ENSG_gene_id)

                gtex_link = 'https://gtexportal.org/home/gene/' + ENSG_gene_id
                omim_link = 'https://www.omim.org/entry/' + omim_id
                gnomAD_link = 'https://gnomad.broadinstitute.org/gene/' + ENSG_gene_id
                decipher_link = 'https://www.deciphergenomics.org/sequence-variant/' + hg38_variant
                clinvar_link = 'https://www.ncbi.nlm.nih.gov/clinvar/?term=' + MANE_select_NM_variant

                uniprot_link, domain_info = get_uniprot_info(ENSG_gene_id)

                expression_eye, expression_brain, expression_fibroblasts, expression_tibial_nerve, expression_blood, expression_transformed_lymphocytes = get_gene_expression(ENSG_gene_id, MANE_select_ENST_variant)

                # to do
                MANE_select_ENST_exon = 'to_do'
                # #partial_lovd_match1 = 'to_do'
                # partial_lovd_match2 = 'to_do'
                # partial_lovd_match3 = 'to_do'
                # partial_lovd_match4 = 'to_do'
                # partial_lovd_match5 = 'to_do'
                # partial_lovd_match6 = 'to_do'
                # partial_lovd_match7 = 'to_do'
                # partial_lovd_match8 = 'to_do'
                # partial_lovd_match9 = 'to_do'
                # partial_lovd_match10 = 'to_do'

                db = get_db()
                db.execute(
                    "INSERT INTO post "
                    # user
                    "(title,"
                    "author_id,"

                    # gene
                    "gene_symbol,"
                    "ENSG_gene_id,"

                    # variant
                    "NC_variant,"
                    "strand,"
                    "hg38_variant,"
                    "MANE_select_NM_variant,"
                    "MANE_select_ENST_variant,"
                    "consequence_variant,"

                    # exon
                    "exon_number,"
                    "total_exons,"
                    "NC_exon,"
                    "exon_length,"
                    "total_protein_length,"
                    "percentage_length,"
                    "frame,"
                    "MANE_select_NM_exon,"
                    "MANE_select_ENST_exon,"
                    "consequence_skipping,"

                    # domain
                    "uniprot_link,"
                    "domain_info,"

                    # expression
                    "gtex_link,"
                    "expression_eye,"
                    "expression_brain,"
                    "expression_fibroblasts,"
                    "expression_tibial_nerve,"
                    "expression_blood,"
                    "expression_transformed_lymphocytes,"

                    # lovd
                    "no_exact_lovd_matches,"
                    "no_partial_lovd_matches,"
                    "exact_lovd_match_link,"
                    "lovd_link,"
                    "partial_lovd_match1,"
                    "partial_lovd_match2,"
                    "partial_lovd_match3,"
                    "partial_lovd_match4,"
                    "partial_lovd_match5,"
                    "partial_lovd_match6,"
                    "partial_lovd_match7,"
                    "partial_lovd_match8,"
                    "partial_lovd_match9,"
                    "partial_lovd_match10,"

                    # identifiers
                    "omim_id,"

                    # links
                    "omim_link,"
                    "gnomAD_link,"
                    "decipher_link,"
                    "clinvar_link)"
                    "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, "
                    "?, ?, ?, ?, ?, ?, ?, ?, ?, ?, "
                    "?, ?, ?, ?, ?, ?, ?, ?, ?, ?,"
                    "?, ?, ?, ?, ?, ?, ?, ?, ?, ?,"
                    "?, ?, ?, ?, ?, ?, ?, ?)",

                    (# user
                    title,
                    g.user["id"],

                    # gene
                    gene_symbol,
                    ENSG_gene_id,

                    # variant
                    NC_variant,
                    strand,
                    hg38_variant,
                    MANE_select_NM_variant,
                    MANE_select_ENST_variant,
                    consequence_variant,

                    # exon
                    exon_number,
                    total_exons,
                    NC_exon,
                    exon_length,
                    total_protein_length,
                    percentage_length,
                    frame,
                    MANE_select_NM_exon,
                    MANE_select_ENST_exon,
                    consequence_skipping,

                    # domain
                    uniprot_link,
                    domain_info,

                    # expression
                    gtex_link,
                    expression_eye,
                    expression_brain,
                    expression_fibroblasts,
                    expression_tibial_nerve,
                    expression_blood,
                    expression_transformed_lymphocytes,

                    # lovd
                    no_exact_lovd_matches,
                    no_partial_lovd_matches,
                    exact_lovd_match_link,
                    lovd_link,
                    partial_lovd_match1,
                    partial_lovd_match2,
                    partial_lovd_match3,
                    partial_lovd_match4,
                    partial_lovd_match5,
                    partial_lovd_match6,
                    partial_lovd_match7,
                    partial_lovd_match8,
                    partial_lovd_match9,
                    partial_lovd_match10,

                    # identifiers
                    omim_id,

                    # links
                    omim_link,
                    gnomAD_link,
                    decipher_link,
                    clinvar_link))

                db.commit()
                return redirect(url_for("blog.output"))
    return render_template("blog/create.html")
