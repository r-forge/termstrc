library("RUnit")

testsuite <- defineTestSuite("termstrc tests",dirs = "/Users/rob/WORK/termstrc/termstrc2/code", testFileRegexp = "unittests.R")

testResult <- runTestSuite(testsuite)

printTextProtocol(testResult)
printHTMLProtocol(testResult, "testResults.html")


