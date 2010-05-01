from selenium import selenium
import unittest, time, re

class SearchPubmed(unittest.TestCase):
    def setUp(self):
        self.verificationErrors = []
        self.selenium = selenium("localhost", 4444, "*chrome", "http://www.ncbi.nlm.nih.gov/pubmed/")
        self.selenium.start()
    
    def test_untitled(self):
        sel = self.selenium
        sel.open("/sites/entrez")
        sel.click("//div[@id='maincontent']/div[4]/div[1]/div[2]/p[3]/span[2]")
        try: self.failUnless(sel.is_text_present("PMID: 20194237"))
        except AssertionError, e: self.verificationErrors.append(str(e))
    
    def tearDown(self):
        self.selenium.stop()
        self.assertEqual([], self.verificationErrors)

if __name__ == "__main__":
    unittest.main()
