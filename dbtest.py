from mysql.connector import connection

# make sure mysql is running on port 3306 - default port
cnx = connection.MySQLConnection(user='root', password='root',
                                 host='localhost',
                                 database='cwg')
cursor = cnx.cursor()
cursor.execute('SELECT * FROM pet')
for row in cursor:
    print row

cursor.close()
cnx.close()
