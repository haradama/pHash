<!DOCTYPE html>
<html>
  <head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">
    <title>pHash Report</title>
    <link rel="stylesheet" type="text/css" href="assets/bootstrap.min.css">
    <style>
      body {
        font-size: 12px;
      }

      .bg-phash {
        background-color: #c8c8c8;
      }
    </style>
  </head>
  <body>
    <section class="bs-docs-section">
      <div class="p-3 mb-2 bg-phash">
        <div class="row">
          <div class="col-md-10 offset-md-1">
            <div class="jumbotron bg-white">
              <div class="col-md-8 offset-md-1">
                <img src="assets/pHash_logo.svg" width="45%" alt="pHash">
              </div>
              <div class="col-xs-12" style="height:30px;"></div>
              <div class="bs-component">
                <table class="table table-hover">
                  <thead>
                    <tr>
                      <th scope="col">Contig ID 
                      </th>
                      <th scope="col">Sequence 
                      </th>
                      <th scope="col">Length 
                      </th>
                      <th scope="col">Similar plasmid 
                      </th>
                      <th scope="col">Phylum 
                      </th>
                      <th scope="col">Jaccard index 
                      </th>
                    </tr>
                  </thead>
                  <tbody>{{range .}} 
                    <tr>
                      <td>{{ .AccID }} 
                      </th>
                      <td>{{ .Seq }} 
                      </th>
                      <td>{{ .Length }} 
                      </th>
                      <td>{{ .Link }} 
                      </th>
                      <td>{{ .Phylum }} 
                      </th>
                      <td>{{ .Jaccard }} 
                      </th>
                    </tr>{{end}} 
                  </tbody>
                </table>
              </div>
            </div>
          </div>
        </div>
      </div>
    </section>
  </body>
  <footer class="footer mt-auto py-3">
    <div class="container">
      <span class="text-muted">https://github.com/haradama/pHash</span>
    </div>
  </footer>
  <script src="assets/bootstrap.min.js"></script>
</html>