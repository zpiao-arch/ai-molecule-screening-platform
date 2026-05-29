import unittest

from fastapi.testclient import TestClient

import server


class StaticAssetsApiTests(unittest.TestCase):
    def setUp(self):
        self.client = TestClient(server.app)

    def test_favicon_request_is_empty_success(self):
        response = self.client.get("/favicon.ico")

        self.assertEqual(response.status_code, 204)
        self.assertEqual(response.content, b"")


if __name__ == "__main__":
    unittest.main()
