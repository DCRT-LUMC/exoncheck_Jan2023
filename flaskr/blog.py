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

bp = Blueprint("blog", __name__)

@bp.route("/")
def welcome():
    """Show the welcome page as default/root """
    return render_template("blog/welcome.html")


@bp.route("/output")
@login_required
def output():
    """ Show all uploaded variants of the current user, most recent first"""
    db = get_db()

    posts = db.execute(
        "SELECT p.id, title, created, author_id, username, "
        # Show gene info
        "gene_symbol,"
        "ENSG_gene_id,"

        # Show variant info
        "NC_variant,"
        "strand,"
        "hg38_variant,"
        "MANE_select_NM_variant,"
        "MANE_select_ENST_variant,"
        "consequence_variant,"

        # Show exon info
        "exon_number,"
        "total_exons,"
        "exon_number_interpretation,"
        "NC_exon,"
        "exon_length,"
        "nearest_splice_distant,"
        "total_protein_length,"
        "percentage_length,"
        "frame,"
        "MANE_select_NM_exon,"
        "MANE_select_ENST_exon,"
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
        " FROM post p JOIN user u ON p.author_id = u.id"
        " WHERE u.id = ? ORDER BY created DESC", (g.user["id"],),
    ).fetchall()
    return render_template("blog/output.html", posts=posts)


def get_post(id, check_author=True):
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
    """Upload new variant by the current user."""
    if request.method == "POST":
        title = request.form["title"]
        syntax_message = check_for_hgvs_format(title)

        # TODO: accept more input formats
        if not title or not title.startswith('NM'):
            flash("Please upload your variant in HGVS format")

        if syntax_message != '':
            flash(syntax_message)

        else:
            # Check if input is in correct HGVS format using VariantValidator, syntax only
            # Check if nucleotides and positions of input correspond to the submitted transcript
            match_message = check_for_match_variant_and_transcript(title)

            if match_message != '':
                flash(match_message)

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
            exon_number_interpretation, \
            NC_exon, \
            exon_length, \
            nearest_splice_distant, \
            total_protein_length, \
            percentage_length, \
            frame, \
            consequence_skipping, \
            MANE_select_NM_exon = exploit_variant_validator(MANE_select_NM_variant)

            # Get LOVD information
            lovd_output = get_lovd_info(hg38_variant, NC_variant)

            # Get strand
            strand = get_strand(ENSG_gene_id)

            # Get links
            gtex_link = 'https://gtexportal.org/home/gene/' + ENSG_gene_id
            omim_link = 'https://www.omim.org/entry/' + omim_id
            gnomAD_link = 'https://gnomad.broadinstitute.org/gene/' + ENSG_gene_id
            decipher_link = 'https://www.deciphergenomics.org/sequence-variant/' + hg38_variant
            clinvar_link = 'https://www.ncbi.nlm.nih.gov/clinvar/?term=' + MANE_select_NM_variant

            uniprot_link, domain_info = get_uniprot_info(ENSG_gene_id)

            # Get expressions
            expression_brain, expression_fibroblasts, expression_tibial_nerve, expression_blood, expression_transformed_lymphocytes = get_gene_expression(ENSG_gene_id, MANE_select_ENST_variant)
            expression_periphery_retina, expression_center_retina = get_eye_expression(ENSG_gene_id)

            # to do
            MANE_select_ENST_exon = 'to_do'


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
                "exon_number_interpretation,"
                "NC_exon,"
                "exon_length,"
                "nearest_splice_distant,"
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
                "VALUES (?, ?, ?, ?, ?, ?, ?,"
                "?, ?, ?, ?, ?, ?, ?, ?, ?, ?, "
                "?, ?, ?, ?, ?, ?, ?, ?, ?, ?,"
                "?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",

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
                exon_number_interpretation,
                NC_exon,
                exon_length,
                nearest_splice_distant,
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
            return redirect(url_for("blog.output"))
    return render_template("blog/create.html")


@bp.route("/<int:id>/update", methods=("GET", "POST"))
@login_required
def update(id):
    """Update a post if the current user is the author."""
    post = get_post(id)

    if request.method == "POST":
        title = request.form["title"]
        body = request.form["body"]
        error = None

        if not title:
            error = "Title is required."

        if error is not None:
            flash(error)
        else:
            db = get_db()
            db.execute(
                "UPDATE post SET title = ?, body = ? WHERE id = ?", (title, body, id)
            )
            db.commit()
            return redirect(url_for("blog.output"))

    return render_template("blog/update.html", post=post)


@bp.route("/<int:id>/delete", methods=("POST",))
@login_required
def delete(id):
    """Delete a post.
    Ensures that the post exists and that the logged in user is the
    author of the post.
    """
    get_post(id)
    db = get_db()
    db.execute("DELETE FROM post WHERE id = ?", (id,))
    db.commit()
    return redirect(url_for("blog.output"))