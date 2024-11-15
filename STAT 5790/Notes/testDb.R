library(DBI)
# dbConnect()
con <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = "mydatabase.db")

DBI::dbListTables(con)

DBI::dbListFields(con, "ecom")

DBI::dbGetQuery(conn = con, statement = "PRAGMA table_info(ecom)")

DBI::dbReadTable(conn = con, name = 'COMPANY')

DBI::dbGetQuery(conn = con, statement = "select * from ecom limit 20")

DBI::dbDisconnect(conn = con)
