library("RUnit")
source("termstrcPackage.R")

testsuite <- defineTestSuite("termstrc tests",dirs = "/home/rferstl/termstrc2/code", testFileRegexp = "unittests.R")

testResult <- runTestSuite(testsuite)

printTextProtocol(testResult)
printHTMLProtocol(testResult, "testResults.html")


