-- Create root password
ALTER USER 'root'@'localhost' IDENTIFIED BY 'rootpass';

-- Allow external connections
DROP USER IF EXISTS 'root'@'%';
CREATE USER 'root'@'%' IDENTIFIED BY 'rootpass';
GRANT ALL ON *.* TO 'root'@'%';
FLUSH PRIVILEGES;
