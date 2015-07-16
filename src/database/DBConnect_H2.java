package database;

import java.io.File;
import java.sql.DriverManager;

public class DBConnect_H2 extends DBConnect {


	/**
	 * 
	 * @param dbPath
	 * @throws Exception
	 */
	public DBConnect_H2(String dbPath) throws Exception {
		connect(dbPath, "sa", "");
	}
	
	
	/**
	 * 
	 * @param dbPath
	 * @param userName
	 * @param password
	 * @throws Exception
	 */
	public DBConnect_H2(String dbPath, String userName, String password) throws Exception {
		connect(dbPath, userName, password);
	}
	
	
	/**
	 * Obtain a connection to the database
	 * @param dbPath
	 * @param userName
	 * @param password
	 * @throws Exception
	 */
	protected void connect(String dbPath, String userName, String password) throws Exception {
		Class.forName("org.h2.Driver");
		this.setConnection(DriverManager.getConnection("jdbc:h2:"+(new File(dbPath)).getAbsolutePath(), userName, password));
		
		this.setQuery_ContainsTable("SELECT table_name FROM INFORMATION_SCHEMA.TABLES where table_schema = 'PUBLIC'");
	}
	
	
	
	
	public static void main(String[] a) throws Exception {
		DBConnect_H2 db = new DBConnect_H2("~/Desktop/TESTDB");
//		db.getStatement().executeQuery("SELECT...");
		db.closeConnection();
	}

}
