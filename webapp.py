import os
import os.path
import sqlite3
import time
import predictor
import io
import uuid
import collections


import cherrypy

from cherrypy.lib.static import serve_fileobj

from jinja2 import Environment, FileSystemLoader

env = Environment(loader=FileSystemLoader('templates'))

DB_STRING = "my.db"


class Root(object):

    @cherrypy.expose
    def index(self):
        tmpl = env.get_template('index.html')

        # allow user to select which algorithm to use

        predictors = predictor.AVAILABLE_PREDICTORS.keys()

        default_pred = 'naive'

        return tmpl.render(predictors = predictors, dp = default_pred)


class GeneViewer(object):

    exposed = True

    def GET(self, uid=None):
        print(uid)
        if uuid is None:
            raise cherrypy.HTTPRedirect("/")

        with sqlite3.connect(DB_STRING) as db:

            tmpl = env.get_template('view.html')

            c = db.cursor()
            c.execute(
                "select scaffold_name, count(*) " +
                "from uuid_gene where uid = ? group by scaffold_name ", (uid,))

            db_genes = c.fetchall()

            # print(genes)

            return tmpl.render(genes=db_genes, uid=uid)


class GFFService(object):

    exposed = True

    def GET(self, uid=None):
        print(uid)
        if uuid is None:
            raise cherrypy.HTTPRedirect("/")

        data = collections.OrderedDict()

        with sqlite3.connect(DB_STRING) as db:

            c = db.cursor()
            c.execute(
                "select scaffold_name, sequence_data " +
                "from uuid_scaffold where uid = ? order by scaffold_name ", (uid,))

            scaffolds = c.fetchall()


            for scaffold in scaffolds:

                data[scaffold[0]] = predictor.ScaffoldData(scaffold[1], [])

                c.execute(
                    "select * from uuid_gene where uid = ? and scaffold_name = ? order by start_pos",
                    (uid, scaffold[0]))

                genes = c.fetchall()

                for gene in genes:                    
                    data[scaffold[0]][1].append(predictor.GeneData(
                        0, gene[3], gene[4], (gene[5], gene[6], gene[7], gene[8]), 
                        (gene[9], gene[10]), gene[11]))

            # print(genes)

            
            c.execute("SELECT tot_gc FROM uuid_tot where uid = ?", (uid,))
            gc_percentage = c.fetchone()

            formatter = predictor.GFF3Formatter()

            gff3_formatted = formatter.format_output(data, gc_percentage[0])

            

            gff3_file = io.StringIO("")


            for line in gff3_formatted:
                print(line, file=gff3_file)

            gff3_file.seek(0)

            return serve_fileobj(
                gff3_file, content_type="text/plain; charset=utf-8", disposition="attachment", name="predictions.gff3")


class StatsViewer(object):

    exposed = True

    def GET(self, uid=None):
        if uid is None:
            raise cherrypy.HTTPRedirect("/")

        with sqlite3.connect(DB_STRING) as db:

            c = db.cursor()
            display_data = []
            gc_percent = 0
            for i in range(4):

                freqs = []
                stmt = "SELECT cast(gc_" + str(i) + "*100 AS INTEGER)  as range_id, count(*) as " + \
                    "range_count from uuid_gene where uid = ? group by range_id"

                c.execute(stmt, (uid,))
                freqs.append(c.fetchall())

                for item in freqs:
                    labels = [x[0] for x in item]
                    values = [x[1] for x in item]
                    display_data.append((labels, values))

                c.execute("SELECT tot_gc FROM uuid_tot where uid = ?", (uid,))
                gc_percent = c.fetchone()

            tmpl = env.get_template('stats.html')

            return tmpl.render(dd=display_data, gc=int(
                gc_percent[0] * 100), uid=uid)


def insert_gene_data(uid, data):

    # handle the data so that it matches what our db wants.

    uid_str = str(uid)

    (gc_pct, scaffolds_data) = data

    db = sqlite3.connect(DB_STRING)

    c = db.cursor()

    c.execute("INSERT INTO uuid_tot VALUES (?, ?)", (uid_str, gc_pct))

    for scaffold_name, contents in scaffolds_data.items():

        c.execute("INSERT INTO uuid_scaffold VALUES(?, ?, ?)",
                  (uid_str, scaffold_name, contents.scaffold_sequence))

        db_format_array = []

        for gene in contents.genes_data:

            (count, start_pos,
             end_pos, gc_content, shine_box, strand) = gene

            db_format_array.append((uid_str, scaffold_name + "_" + str(count) + "_" + str(strand),
                                    scaffold_name,
                                    start_pos, end_pos,
                                    *gc_content, *shine_box, strand))

        c.executemany(
            "INSERT INTO uuid_gene VALUES (?, ?, ?, ?, ?, ?, ?, ?, ? ,? ,?, ?)",
            db_format_array)

    db.commit()


class DNAUploadService(object):
    exposed = True

    def POST(self, myFile, selected_predictor, shortest_gene, intra_gene_gap,
             shine_box_distance, **kwargs):

        allow_runoff = False

        if 'runoff' in kwargs:
            allow_runoff = True

        uid = uuid.uuid4()
        with sqlite3.connect(DB_STRING) as c:
            cherrypy.session['ts'] = time.time()

            c.execute("INSERT INTO user_uuid VALUES (?, ?)",
                      [cherrypy.session.id, str(uid)])

        config = predictor.GpConfig(
            selected_predictor,
            int(shortest_gene),
            allow_runoff,
            int(intra_gene_gap),
            int(shine_box_distance))

        try:
            gene_data = predictor.handleInput(
                io.TextIOWrapper(myFile.file), config)
        except:
            raise cherrypy.HTTPRedirect("/")

        insert_gene_data(uid, gene_data)

        raise cherrypy.HTTPRedirect("/stats/" + str(uid))


def setup_database():
    """
    Create the `user_string` table in the database
    on server startup
    """
    with sqlite3.connect(DB_STRING) as con:
        con.execute("CREATE TABLE user_uuid (session_id, uid)")
        con.execute("CREATE TABLE uuid_tot (uid TEXT, tot_gc REAL)")
        con.execute("CREATE TABLE uuid_scaffold (uid, scaffold_name TEXT,  sequence_data TEXT, \
                                          PRIMARY KEY (uid, scaffold_name))")
        con.execute("CREATE TABLE uuid_gene (uid, gene_name TEXT, scaffold_name TEXT, \
                                          start_pos INTEGER, end_pos INTEGER, gc_0 REAL, gc_1 REAL, gc_2 REAL, gc_3 REAL, \
                                          shine_str TEXT, shine_pos INTEGER, strand INTEGER, PRIMARY KEY (uid, gene_name), \
                                          FOREIGN KEY (uid, scaffold_name) REFERENCES uuid_scaffold(uid, scaffold_name))")


def cleanup_database():
    """
    Destroy the `user_string` table from the database
    on server shutdown.
    """
    with sqlite3.connect(DB_STRING) as con:
        con.execute("DROP TABLE user_uuid")
        con.execute("DROP TABLE uuid_gene")
        con.execute("DROP TABLE uuid_scaffold")
        con.execute("DROP TABLE uuid_tot")

if __name__ == '__main__':
    conf = {
        '/': {
            'tools.sessions.on': True,
        },
        '/dna': {
            'request.dispatch': cherrypy.dispatch.MethodDispatcher(),
            'tools.sessions.on': True,
            'tools.response_headers.on': True,
            'tools.response_headers.headers': [('Content-Type', 'text/plain')],
        },
        '/gff3': {
            'request.dispatch': cherrypy.dispatch.MethodDispatcher(),
            'tools.sessions.on': True,
            'tools.response_headers.on': True,
        },
        '/view': {
            'tools.sessions.on': True,
            'request.dispatch': cherrypy.dispatch.MethodDispatcher(),
        },
        '/stats': {
            'tools.sessions.on': True,
            'request.dispatch': cherrypy.dispatch.MethodDispatcher(),
        },

    }

    cherrypy.engine.subscribe('start', setup_database)
    cherrypy.engine.subscribe('stop', cleanup_database)

    webapp = Root()
    webapp.dna = DNAUploadService()
    webapp.gff3 = GFFService()
    webapp.view = GeneViewer()
    webapp.stats = StatsViewer()
    cherrypy.server.socket_host = '0.0.0.0'
    cherrypy.quickstart(webapp, '/', conf)
