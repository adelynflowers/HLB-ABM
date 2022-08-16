#include <stdlib.h>
#include <iostream>
#include <mysql_connection.h>
#include <driver.h>
#include <exception.h>
#include <resultset.h>
#include <statement.h>

using namespace sql;
using namespace std;
int main(void){
  sql::Driver *driver;
  sql::Connection *con;
  sql::Statement *stmt;
  sql::ResultSet *res;

  driver = get_driver_instance();
  con = driver->connect("tcp://127.0.0.1:3306","sam","CitrusABM21");
  con->setSchema("citrus");
  stmt = con->createStatement();
  res = stmt->executeQuery("SELECT 'Hello World!' AS _message");
  while (res->next()) {
      cout << "My SQL replies: ";
      cout << res->getString("_message") << endl;
      cout << "MYSQL says it again: ";
      cout << res->getString(1) << endl;
  }
  delete res;
  delete stmt;
  delete con; 
  return 0;
}